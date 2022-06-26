
## The purpose of this script is to implement the pentanucleotide model for the pathway dN/dS in psoriasis. 

.libPaths("/lustre/scratch119/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")
library("MASS")
library("GenomicRanges")
library("Rsamtools")
library("seqinr")
library(dndscv)

load("/nfs/team78pc20/im3/Projects/Selection_paper/Pentanucleotide_models/Lmatrix_pentanucleotide_pergene_withoutnormalpanelbias.RData")
genomeFile = "/nfs/cancer_ref01/Homo_sapiens/37/genome.fa"
output_directory = "/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/pathway_pentamodel/"

## I used the USCS liftOver tool to lift over the mutations. They were originally called on hg38. 
mutations_file="/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/all_mutations_hg19.txt"
#passengers = read.table("/nfs/team78pc16/im3/Reference_data/putative_passengers_from_multiple_sources_20170606.txt", header=0, sep="\t", stringsAsFactors=F)[,1]
passengers = read.table("/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/passenger_genes_pathway_dNdS.txt", header=0, sep="\t", stringsAsFactors=F)[,1]
g = read.table("/nfs/team78pc16/im3/Reference_data/Reference_list_genes_Human_CDS_dNdS_new.txt", header=0, sep="\t", stringsAsFactors=F)[,1]

parameters <- commandArgs(TRUE)
geneList  <- parameters[1]
geneL <- read.table(paste("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/08_selection_analyses/pathway_dnds_gene_lists/",geneList, ".txt", sep=""))
drivers <- geneL$V1

## Define variables that will be the same for all gene sets. Only needs doing once. 
mutations = read.table(mutations_file, header=1, sep="\t", stringsAsFactors=F)
nt = c("A","C","G","T")
compnt = rev(nt); names(compnt) = nt

L_matrix_bck <- L_matrix

annotate_mutations <- function (mutations, gene_list = NULL, refdb = "hg19", sm = "192r_3w", 
                                kc = "cgc81", cv = "hg19", max_muts_per_gene_per_sample = Inf, 
                                max_coding_muts_per_sample = Inf, use_indel_sites = T, min_indels = 5, 
                                maxcovs = 20, constrain_wnon_wspl = T, outp = 3) 
{
  message("[1] Loading the environment...")
  cat("max_muts_per_gene_per_sample=",max_muts_per_gene_per_sample,"\n")
  cat("max_coding_muts_per_sample=",max_coding_muts_per_sample,"\n")
  if (refdb == "hg19") {
    data("refcds_hg19", package = "dndscv")
  }
  else {
    load(refdb)
  }
  if (is.null(gene_list)) {
    gene_list = sapply(RefCDS, function(x) x$gene_name)
  }
  else {
    allg = sapply(RefCDS, function(x) x$gene_name)
    nonex = gene_list[!(gene_list %in% allg)]
    if (length(nonex) > 0) {
      stop(sprintf("The following input gene names are not in the RefCDS database: %s", 
                   paste(nonex, collapse = ", ")))
    }
    RefCDS = RefCDS[allg %in% gene_list]
    gr_genes = gr_genes[gr_genes$names %in% gene_list]
  }
  if (is.character(cv)) {
    data(list = sprintf("covariates_%s", cv), package = "dndscv")
  }
  else {
    covs = cv
  }
  if (kc[1] %in% c("cgc81")) {
    data(list = sprintf("cancergenes_%s", kc), package = "dndscv")
  }
  else {
    known_cancergenes = kc
  }
  if (length(sm) == 1) {
    data(list = sprintf("submod_%s", sm), package = "dndscv")
  }
  else {
    substmodel = sm
  }
  for (j in 1:length(RefCDS)) {
    RefCDS[[j]]$seq_cds = base::strsplit(as.character(RefCDS[[j]]$seq_cds), 
                                         split = "")[[1]]
    RefCDS[[j]]$seq_cds1up = base::strsplit(as.character(RefCDS[[j]]$seq_cds1up), 
                                            split = "")[[1]]
    RefCDS[[j]]$seq_cds1down = base::strsplit(as.character(RefCDS[[j]]$seq_cds1down), 
                                              split = "")[[1]]
    if (!is.null(RefCDS[[j]]$seq_splice)) {
      RefCDS[[j]]$seq_splice = base::strsplit(as.character(RefCDS[[j]]$seq_splice), 
                                              split = "")[[1]]
      RefCDS[[j]]$seq_splice1up = base::strsplit(as.character(RefCDS[[j]]$seq_splice1up), 
                                                 split = "")[[1]]
      RefCDS[[j]]$seq_splice1down = base::strsplit(as.character(RefCDS[[j]]$seq_splice1down), 
                                                   split = "")[[1]]
    }
  }
  message("[2] Annotating the mutations...")
  colnames(mutations) = c("sampleID", "chr", "pos", "ref", 
                          "mut")
  nt = c("A", "C", "G", "T")
  trinucs = paste(rep(nt, each = 16, times = 1), rep(nt, each = 4, 
                                                     times = 4), rep(nt, each = 1, times = 16), sep = "")
  trinucinds = setNames(1:64, trinucs)
  trinucsubs = NULL
  for (j in 1:length(trinucs)) {
    trinucsubs = c(trinucsubs, paste(trinucs[j], paste(substr(trinucs[j], 
                                                              1, 1), setdiff(nt, substr(trinucs[j], 2, 2)), substr(trinucs[j], 
                                                                                                                   3, 3), sep = ""), sep = ">"))
  }
  trinucsubsind = setNames(1:192, trinucsubs)
  ind = setNames(1:length(RefCDS), sapply(RefCDS, function(x) x$gene_name))
  gr_genes_ind = ind[gr_genes$names]
  if (any(diff(mutations$pos) == 1)) {
    warning("Mutations observed in contiguous sites within a sample. Please annotate or remove dinucleotide or complex substitutions for best results.")
  }
  if (nrow(unique(mutations[, 2:5])) < nrow(mutations)) {
    warning("Same mutations observed in different sampleIDs. Please verify that these are independent events and remove duplicates otherwise.")
  }
  gr_muts = GRanges(mutations$chr, IRanges(mutations$pos, mutations$pos))
  ol = as.matrix(findOverlaps(gr_muts, gr_genes, type = "any", 
                              select = "all"))
  mutations = mutations[ol[, 1], ]
  mutations$geneind = gr_genes_ind[ol[, 2]]
  mutations$gene = sapply(RefCDS, function(x) x$gene_name)[mutations$geneind]
  nsampl = sort(table(mutations$sampleID))
  exclsamples = NULL
  if (any(nsampl > max_coding_muts_per_sample)) {
    message(sprintf("    Note: %0.0f samples excluded for exceeding the limit of mutations per sample", 
                    sum(nsampl > max_coding_muts_per_sample)))
    exclsamples = names(nsampl[nsampl > max_coding_muts_per_sample])
    mutations = mutations[!(mutations$sampleID %in% names(nsampl[nsampl > 
                                                                   max_coding_muts_per_sample])), ]
  }
  mutrank = ave(mutations$pos, paste(mutations$sampleID, mutations$gene), 
                FUN = function(x) rank(x))
  exclmuts = NULL
  if (any(mutrank > max_muts_per_gene_per_sample)) {
    message(sprintf("    Note: %0.0f mutations removed for exceeding the limit of mutations per gene per sample", 
                    sum(mutrank > max_muts_per_gene_per_sample)))
    exclmuts = mutations[mutrank > max_muts_per_gene_per_sample, 
    ]
    mutations = mutations[mutrank <= max_muts_per_gene_per_sample, 
    ]
  }
  snv = (mutations$ref %in% nt & mutations$mut %in% nt)
  indels = mutations[!snv, ]
  mutations = mutations[snv, ]
  mutations$ref_cod = mutations$ref
  mutations$mut_cod = mutations$mut
  compnt = setNames(rev(nt), nt)
  mutations$strand = sapply(RefCDS, function(x) x$strand)[mutations$geneind]
  isminus = (mutations$strand == -1)
  mutations$ref_cod[isminus] = compnt[mutations$ref[isminus]]
  mutations$mut_cod[isminus] = compnt[mutations$mut[isminus]]
  for (j in 1:length(RefCDS)) {
    RefCDS[[j]]$N = array(0, dim = c(192, 4))
  }
  chr2cds = function(pos, cds_int, strand) {
    if (strand == 1) {
      return(which(pos == unlist(apply(cds_int, 1, function(x) x[1]:x[2]))))
    }
    else if (strand == -1) {
      return(which(pos == rev(unlist(apply(cds_int, 1, 
                                           function(x) x[1]:x[2])))))
    }
  }
  ref3_cod = mut3_cod = wrong_ref = aachange = ntchange = impact = array(NA, 
                                                                         nrow(mutations))
  for (j in 1:nrow(mutations)) {
    geneind = mutations$geneind[j]
    pos = mutations$pos[j]
    if (any(pos == RefCDS[[geneind]]$intervals_splice)) {
      impact[j] = "Essential_Splice"
      impind = 4
      pos_ind = (pos == RefCDS[[geneind]]$intervals_splice)
      cdsnt = RefCDS[[geneind]]$seq_splice[pos_ind]
      ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], 
                            RefCDS[[geneind]]$seq_splice[pos_ind], RefCDS[[geneind]]$seq_splice1down[pos_ind])
      mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], 
                            mutations$mut_cod[j], RefCDS[[geneind]]$seq_splice1down[pos_ind])
      aachange[j] = ntchange[j] = "-"
    }
    else {
      pos_ind = chr2cds(pos, RefCDS[[geneind]]$intervals_cds, 
                        RefCDS[[geneind]]$strand)
      cdsnt = RefCDS[[geneind]]$seq_cds[pos_ind]
      ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], 
                            RefCDS[[geneind]]$seq_cds[pos_ind], RefCDS[[geneind]]$seq_cds1down[pos_ind])
      mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], 
                            mutations$mut_cod[j], RefCDS[[geneind]]$seq_cds1down[pos_ind])
      codon_pos = c(ceiling(pos_ind/3) * 3 - 2, ceiling(pos_ind/3) * 
                      3 - 1, ceiling(pos_ind/3) * 3)
      old_codon = as.character(as.vector(RefCDS[[geneind]]$seq_cds[codon_pos]))
      pos_in_codon = pos_ind - (ceiling(pos_ind/3) - 1) * 
        3
      new_codon = old_codon
      new_codon[pos_in_codon] = mutations$mut_cod[j]
      old_aa = seqinr::translate(old_codon)
      new_aa = seqinr::translate(new_codon)
      aachange[j] = sprintf("%s%0.0f%s", old_aa, ceiling(pos_ind/3), 
                            new_aa)
      ntchange[j] = sprintf("%s%0.0f%s", mutations$ref_cod[j], 
                            pos_ind, mutations$mut_cod[j])
      if (new_aa == old_aa) {
        impact[j] = "Synonymous"
        impind = 1
      }
      else if (new_aa == "*") {
        impact[j] = "Nonsense"
        impind = 3
      }
      else if (old_aa != "*") {
        impact[j] = "Missense"
        impind = 2
      }
      else if (old_aa == "*") {
        impact[j] = "Stop_loss"
        impind = NA
      }
    }
    if (mutations$ref_cod[j] != as.character(cdsnt)) {
      wrong_ref[j] = 1
    }
    else if (!is.na(impind)) {
      trisub = trinucsubsind[paste(ref3_cod[j], mut3_cod[j], 
                                   sep = ">")]
      RefCDS[[geneind]]$N[trisub, impind] = RefCDS[[geneind]]$N[trisub, 
                                                                impind] + 1
    }
    if (round(j/10000) == (j/10000)) {
      message(sprintf("    %0.3g %%...", round(j/nrow(mutations), 
                                               2) * 100))
    }
  }
  mutations$ref3_cod = ref3_cod
  mutations$mut3_cod = mut3_cod
  mutations$aachange = aachange
  mutations$ntchange = ntchange
  mutations$impact = impact
  mutations$pid = sapply(RefCDS, function(x) x$protein_id)[mutations$geneind]
  if (any(!is.na(wrong_ref))) {
    stop(sprintf("%0.0f mutations have a wrong reference base, please correct and rerun.", 
                 sum(!is.na(wrong_ref))))
    wrong_refbase = mutations[!is.na(wrong_ref), ]
    mutations = mutations[is.na(wrong_ref), ]
  }
  if (any(nrow(indels))) {
    indels = cbind(indels, data.frame(ref_cod = ".", mut_cod = ".", 
                                      strand = ".", ref3_cod = ".", mut3_cod = ".", aachange = ".", 
                                      ntchange = ".", impact = "no-SNV", pid = sapply(RefCDS, 
                                                                                      function(x) x$protein_id)[indels$geneind]))
    annot = rbind(mutations, indels)
  }
  else {
    annot = mutations
  }
  annot = annot[order(annot$sampleID, annot$chr, annot$pos), 
  ]
  
  # Por aquí devolver el objeto!
  return(mutations);
}


mutations <- annotate_mutations(mutations,max_muts_per_gene_per_sample=Inf,max_coding_muts_per_sample=Inf)

## Extract the pentamers for each mutation
pentant = as.vector(scanFa(genomeFile,GRanges(mutations$chr, IRanges(mutations$pos-2,mutations$pos+2))))
pentant[mutations$strand==-1] = sapply(pentant[mutations$strand==-1], function(x) paste(compnt[rev(strsplit(x,split="")[[1]])],collapse="")) # Reverse comp for mutations in genes in -1 strand
mutations$ref5_cod = pentant
mutations$mut5_cod = paste(substr(pentant,1,2), mutations$mut_cod, substr(pentant,4,5), sep="")
mutations$penta_sub = paste(mutations$ref5_cod, mutations$mut5_cod, sep=">")


pentanuc_list = paste(rep(nt,each=4^4,times=1),rep(nt,each=4^3,times=4),rep(nt,each=4^2,times=4^2),rep(nt,each=4,times=4^3),rep(nt,each=1,times=4^4), sep="")
pentanuc_subs = c(sapply(pentanuc_list, function(x) paste(x, paste(substr(x,1,2),nt[nt!=substr(x,3,3)],substr(x,4,5),sep=""), sep=">")))
pentanuc_subs_ind = 1:3072; names(pentanuc_subs_ind) = pentanuc_subs

# 3. Creating the N matrix

get_N_matrix = function(m) {
  N_matrix = array(0, dim=c(3072,4))
  impind = 1:4; names(impind) = c("Synonymous","Missense","Nonsense","Essential_Splice")
  imp = impind[m$impact]
  matrind = as.numeric(pentanuc_subs_ind[m$penta_sub])
  # Annotating the counts for synonymous, missense and nonsense mutations
  for (h in 1:4) {
    matrix_ind = table(matrind[which(imp==h)])
    N_matrix[as.numeric(names(matrix_ind)),h] = matrix_ind
  }
  return(N_matrix)
}



Np = get_N_matrix(m = mutations[mutations$gene %in% passengers,])
Lp = apply(L_matrix_bck[,,which(g %in% passengers)], c(1,2), sum)

# 4. MLEs of all parameters using Poisson regression

# Subfunction: Poisson regression given R_matrix, L_matrix and N_matrix

get_mlerates = function(R_matrix, L_matrix, N_matrix, prefixfn, wparams=c("wMIS","wNON","wSPL")) {
  
  R = c(R_matrix); N = c(N_matrix); L = c(L_matrix)
  rmpos = (L==0); R = R[!rmpos]; N = N[!rmpos]; L = L[!rmpos]
  predictors = unique(strsplit(x=paste(R,collapse="*"), split="\\*")[[1]])
  
  indicator_matrix = array(0,dim=c(length(R),length(predictors)))
  predictors_ind = 1:length(predictors); names(predictors_ind) = predictors
  aux = sapply(strsplit(R,"\\*"), function(x) as.numeric(predictors_ind[x]))
  for (j in 1:length(R)) {
    indicator_matrix[j,aux[[j]]] = 1
  }
  indicator_matrix = as.data.frame(indicator_matrix)
  colnames(indicator_matrix) = predictors
  model = glm(formula = N ~ offset(log(L)) + . -1, data=indicator_matrix, family=poisson(link=log))
  write.table(c(logLik=as.vector(logLik(model)),df=df.residual(model)), file=sprintf('%s/%s_logLik.txt',output_directory,prefixfn), row.names = T, col.names = F, sep="\t", quote=F)
  ci = exp(confint.default(model,parm=wparams)) # Wald confidence intervals
  params = exp(coefficients(model))[wparams]
  
  poissreg = list(N=N,L=L,R=R,indicator_matrix=indicator_matrix,model=model,params=params,ci=ci)
  #save(poissreg, file=sprintf("%s/%s_poissreg.RData",output_directory,prefixfn) )
  CItable = data.frame(omega=wparams, MLEs=exp(coefficients(model))[wparams], lowbd=ci[wparams,1], highbd=ci[wparams,2], P=coef(summary(model))[wparams,4])
  write.table(CItable,file=sprintf('%s/%s_dNdSvals.txt',output_directory,prefixfn), row.names = F, col.names = T, sep="\t", quote=F)
  save(model, file=sprintf('%s/%s_glm.rda',output_directory,prefixfn))
  return(CItable)
}



  ## Define variables for drivers
  Ld = apply(L_matrix_bck[,,which(g %in% drivers)], c(1,2), sum)
  Nd = get_N_matrix(m = mutations[mutations$gene %in% drivers,])
 
  L_matrix = rbind(Ld,Lp)
  N_matrix = rbind(Nd,Np)
  
  R_matrix = array("",dim=dim(L_matrix))
  R_matrix[,1] = c(paste(pentanuc_subs,"*r_drivpass",sep="") , pentanuc_subs) # This adds an extra parameter to account for a different mutation rate (synonymous density) for drivers and passengers
  R_matrix[,2] = c(paste(pentanuc_subs,"*r_drivpass*wmis_driv",sep="") , paste(pentanuc_subs,"*wmis_pass",sep=""))
  R_matrix[,3] = c(paste(pentanuc_subs,"*r_drivpass*wnon_driv",sep="") , paste(pentanuc_subs,"*wnon_pass",sep=""))
  R_matrix[,4] = c(paste(pentanuc_subs,"*r_drivpass*wspl_driv",sep="") , paste(pentanuc_subs,"*wspl_pass",sep=""))
  
  dndsout_penta = get_mlerates(R_matrix, L_matrix, N_matrix, paste(geneList, "_Full3075_1x2w_model", sep=""), wparams=c("wmis_driv","wnon_driv","wspl_driv","wmis_pass","wnon_pass","wspl_pass","r_drivpass"))
  write.table(dndsout_penta, file=paste(output_directory, geneList, "_Full3075_1x2w_model.txt", sep=""), sep="\t", quote=F, row.names=F)





