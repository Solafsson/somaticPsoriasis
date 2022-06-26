
.libPaths("/lustre/scratch119/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")

## The following dNdS analyses will be run:
## 1. Unbiased gene-level dN/dS of the entire dataset. (Check the impact of using the hg38 covariates)
## also test the impact of running indels and DBS separately. 
## Test the impact of excluding some outliers (which contribute large number of passengers and comparatively few drivers). 
## 2. Unbiased gene-level dN/dS of mutations found in lesional skin only.
## 3. Pathway-level dN/dS of pathways defined a-priory. 
## 4. Unbiased sitednds, codondnds and domaindnds
## 5. Further restricted hypothesis testing of previously reported genes that don't reach significance in this analysis. 

library("seqinr")
library("Biostrings")
library("MASS")
library("GenomicRanges")
library("dndscv")

sample_meta <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/sample_info/sample_meta.txt", h=T)
patient_list <- unique(sample_meta$patient_ID[sample_meta$Seq_status=="WES" & is.na(sample_meta$excl_criteria)])

patient_list <- c("patient01", "patient03","patient04","patient05", "patient06", "patient07","patient09","patient11",
                  "patient13", "patient17", "patient23", "patient29", "patient31", "patient32","patient34",
                  "patient02","patient10","patient12","patient14","patient16","patient18","patient19",                "patient20","patient21","patient22","patient24","patient26","patient27","patient28","patient30",
                  "patient08","patient115","patient37","patient46","patient55","patient76", "patient97",
                  "patient38","patient47","patient56","patient79","patient98","patient109",
                  "patient29","patient39","patient48","patient57","patient99","patient102","patient40",
                  "patient49","patient58","patient80","patient103","patient110","patient41","patient50",
                  "patient60","patient84","patient104","patient42","patient51","patient62","patient86",
                  "patient112","patient25","patient43","patient52","patient70","patient87","patient105",
                  "patient106","patient113","patient36","patient44","patient53","patient95","patient111")

samples <- as.character()
sample_muts <- as.character()
for(patient in patient_list) {
  #genotype_bin <- read.table(paste("/Users/so11/phd/psoriasis/scratch/binomial_filters/", patient, "/",patient,"_genotype_allMuts.txt", sep=""), h=T)
  genotype_bin <- read.table(paste("/lustre/scratch119/humgen/projects/psoriasis/binomial_filters/", patient, "/",patient,"_genotype_allMuts.txt", sep=""), h=T)
  
  
  ## I haven't made good trees yet so as a temporary solution I am going to keep the mutation
  ## only in the sample in which it has the highest VAF. That way I don't count the same mutation twice. 
  
  for(i in 1:nrow(genotype_bin)) {
    samples <- append(samples, paste(colnames(genotype_bin)[genotype_bin[i,]>0], collapse=","))
  }
  
  sample_muts <- append(sample_muts, rownames(genotype_bin))
  
}

all_mutations <- data.frame(SampleID=samples,
                            Chr=unlist(strsplit(sample_muts, split=":"))[c(T,F,F,F)],
                            Pos=unlist(strsplit(sample_muts, split=":"))[c(F,T,F,F)],
                            Ref=unlist(strsplit(sample_muts, split=":"))[c(F,F,T,F)],
                            Alt=unlist(strsplit(sample_muts, split=":"))[c(F,F,F,T)])
all_mutations <- subset(all_mutations, all_mutations$SampleID!="")


## or
all_mutations <- read.table("/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/all_mutations.txt", h=T)
all_mutations$Chr <- gsub("chr", "", all_mutations$Chr)

## Download the hg38 reference rda file from here:
## https://github.com/im3sanger/dndscv_data/blob/master/data/RefCDS_human_GRCh38.p12.rda

## Main Unbiased analysis
################################
## dNdScv lumps DBS mutations together with indels. This leads to inflated dNdS ratios for indels 
## and can potentially lead to false positive calls.
## There are 54537 DBS mutations in the dataset compared with 3670 indels to the DBS mutations dominate
## the indel signal. (Nearly 15 times more DBS mutations). They have vastly different overdispersion values. 

covs = "/lustre/scratch119/humgen/projects/psoriasis/resources/covariates_20pc_GRCh37-38.altogether_noepiout.Rdat"
load(covs) # it loads an object called scores
refcds_38 = "/lustre/scratch119/humgen/projects/psoriasis/resources/refcds_GRCh38-GencodeV18+Appris.rda"

d38 = dndscv(all_mutations, refdb=refcds_38, cv=scores,max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf) 

sel_cv = d38$sel_cv
signif_genes = sel_cv[sel_cv$qglobal_cv<0.05 | sel_cv$qallsubs_cv<0.05, c(1:6, 17,19)]
signif_genes

## I want to split the "no-SNV" mutation types into indels and DBS, run a separate binomial model
## for each type and then combine them with the P-values from the SNVs. 

dbs_only <- all_mutations[(nchar(all_mutations$Ref)==1 & nchar(all_mutations$Alt)==1) | ( nchar(all_mutations$Ref)==2 & nchar(all_mutations$Alt)==2),]
d38_dbs_only = dndscv(dbs_only, refdb=refcds_38, cv=scores,max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf) 
sel_cv_dbs = d38_dbs_only$sel_cv

indel_only <- all_mutations[nchar(all_mutations$Ref)!=2 | nchar(all_mutations$Alt)!=2,]
d38_indel_only = dndscv(indel_only, refdb=refcds_38, cv=scores,max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf) 
sel_cv_indel = d38_indel_only$sel_cv

# Fisher combined p-values (single base substitutions, double base substitutions and indels)
# sort the dataframes in the same way:
sel_cv_indel = sel_cv_indel[order(sel_cv_indel$gene_name),]
sel_cv_dbs = sel_cv_dbs[order(sel_cv_dbs$gene_name),]

p_global <- 1-pchisq(-2 * (log(sel_cv_indel$pallsubs_cv) + log(sel_cv_indel$pind_cv) + log(sel_cv_dbs$pind_cv)), df = 6)
q_global <- p.adjust(p_global, method="BH")

sel_cv_indel$pdbs_cv <- sel_cv_dbs$pind_cv
sel_cv_indel$wdbs_cv <- sel_cv_dbs$wind_cv
sel_cv_indel$pglobal_cv <- p_global
sel_cv_indel$qglobal_cv <- q_global
sel_cv_indel$n_dbs <- sel_cv_dbs$n_ind

sel_cv_indel = sel_cv_indel[order(sel_cv_indel$pglobal_cv, sel_cv_indel$pallsubs_cv, sel_cv_indel$pmis_cv, sel_cv_indel$ptrunc_cv, -sel_cv_indel$wmis_cv),] # Sorting genes in the output file
                           
mutation_list <- d38$annotmuts
mutations_in_signif_genes <- mutation_list[mutation_list$gene %in% signif_genes$gene_name,]

write.table(mutation_list, "/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/annotated_coding_mutations.txt", sep="\t", row.names=F, quote=F)
write.table(sel_cv_indel, "/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/sel_cv_indel_test.txt", sep="\t", row.names=F, quote=F)

##############################################################


## Test the impact of not using the covariates.  
dbs_noCov = dndscv(dbs_only, refdb=refcds_38, cv=NULL,max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf) 

indel_noCov = dndscv(indel_only, refdb=refcds_38, cv=NULL,max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf) 


sel_cv_dbs_noCov = dbs_noCov$sel_cv
sel_cv_indel_noCov = indel_noCov$sel_cv

sel_cv_indel_noCov = sel_cv_indel_noCov[order(sel_cv_indel_noCov$gene_name),]
sel_cv_dbs_noCov = sel_cv_dbs_noCov[order(sel_cv_dbs_noCov$gene_name),]

p_global <- 1-pchisq(-2 * (log(sel_cv_indel_noCov$pallsubs_cv) + log(sel_cv_indel_noCov$pind_cv) + log(sel_cv_dbs_noCov$pind_cv)), df = 6)
q_global <- p.adjust(p_global, method="BH")

sel_cv_indel_noCov$pdbs_cv <- sel_cv_dbs_noCov$pind_cv
sel_cv_indel_noCov$wdbs_cv <- sel_cv_dbs_noCov$wind_cv
sel_cv_indel_noCov$pglobal_cv <- p_global
sel_cv_indel_noCov$qglobal_cv <- q_global
sel_cv_indel_noCov$n_dbs <- sel_cv_dbs_noCov$n_ind

sel_cv_indel_noCov = sel_cv_indel_noCov[order(sel_cv_indel_noCov$pglobal_cv, sel_cv_indel_noCov$pallsubs_cv, sel_cv_indel_noCov$pmis_cv, sel_cv_indel_noCov$ptrunc_cv, -sel_cv_indel_noCov$wmis_cv),] # Sorting genes in the output file

write.table(sel_cv_indel_noCov, "/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/sel_cv_nocov.txt", sep="\t", row.names=F, quote=F)



## See if the results are impacted by exclusion of hypermutators (which contribute many passengers but few drivers)
## This involves doing two things: First, use the default settings for max_muts_per_gene and sample. Second, actually
## remove hypermutators from the dataset. 

dbs_noHype <- dbs_only[grep("P34H", dbs_only$SampleID, invert=T),]
indel_noHype <- indel_only[grep("P34H", indel_only$SampleID, invert=T),]

dbs_hype = dndscv(dbs_noHype, refdb=refcds_38, cv=scores) 
indel_hype = dndscv(indel_noHype, refdb=refcds_38, cv=scores) 

sel_cv_dbs_hype = dbs_hype$sel_cv
sel_cv_indel_hype = indel_hype$sel_cv

sel_cv_indel_hype = sel_cv_indel_hype[order(sel_cv_indel_hype$gene_name),]
sel_cv_dbs_hype = sel_cv_dbs_hype[order(sel_cv_dbs_hype$gene_name),]

p_global <- 1-pchisq(-2 * (log(sel_cv_indel_hype$pallsubs_cv) + log(sel_cv_indel_hype$pind_cv) + log(sel_cv_dbs_hype$pind_cv)), df = 6)
q_global <- p.adjust(p_global, method="BH")

sel_cv_indel_hype$pdbs_cv <- sel_cv_dbs_hype$pind_cv
sel_cv_indel_hype$wdbs_cv <- sel_cv_dbs_hype$wind_cv
sel_cv_indel_hype$pglobal_cv <- p_global
sel_cv_indel_hype$qglobal_cv <- q_global
sel_cv_indel_hype$n_dbs <- sel_cv_dbs_hype$n_ind

sel_cv_indel_hype = sel_cv_indel_hype[order(sel_cv_indel_hype$pglobal_cv, sel_cv_indel_hype$pallsubs_cv, sel_cv_indel_hype$pmis_cv, sel_cv_indel_hype$ptrunc_cv, -sel_cv_indel_hype$wmis_cv),] # Sorting genes in the output file

write.table(sel_cv_indel_hype, "/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/sel_cv_noHype.txt", sep="\t", row.names=F, quote=F)




## Try running the analysis excluding mutations found in non-lesional skin:
lesional_mutations <- all_mutations[!grepl("H", all_mutations$SampleID),]

dbs_lesional <- lesional_mutations[(nchar(lesional_mutations$Ref)==1 & nchar(lesional_mutations$Alt)==1) | ( nchar(lesional_mutations$Ref)==2 & nchar(lesional_mutations$Alt)==2),]
d38_dbs_lesional = dndscv(dbs_lesional, refdb=refcds_38, cv=scores,max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf) 
sel_cv_dbs = d38_dbs_lesional$sel_cv

indel_lesional <- lesional_mutations[nchar(lesional_mutations$Ref)!=2 | nchar(lesional_mutations$Alt)!=2,]
d38_indel_lesional = dndscv(indel_lesional, refdb=refcds_38, cv=scores,max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf) 
sel_cv_indel = d38_indel_lesional$sel_cv


sel_cv_indel = sel_cv_indel[order(sel_cv_indel$gene_name),]
sel_cv_dbs = sel_cv_dbs[order(sel_cv_dbs$gene_name),]

p_global <- 1-pchisq(-2 * (log(sel_cv_indel$pallsubs_cv) + log(sel_cv_indel$pind_cv) + log(sel_cv_dbs$pind_cv)), df = 6)
q_global <- p.adjust(p_global, method="BH")

sel_cv_indel$pdbs_cv <- sel_cv_dbs$pind_cv
sel_cv_indel$wdbs_cv <- sel_cv_dbs$wind_cv
sel_cv_indel$pglobal_cv <- p_global
sel_cv_indel$qglobal_cv <- q_global
sel_cv_indel$n_dbs <- sel_cv_dbs$n_ind

sel_cv_indel = sel_cv_indel[order(sel_cv_indel$pglobal_cv, sel_cv_indel$pallsubs_cv, sel_cv_indel$pmis_cv, sel_cv_indel$ptrunc_cv, -sel_cv_indel$wmis_cv),] # Sorting genes in the output file
                           
write.table(sel_cv_indel, "/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/sel_cv_lesional.txt", sep="\t", row.names=F, quote=F)



###
##
##      PATHWAY LEVEL dNdS
####################################################

source("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/08_selection_analyses/dndscv_pvalpois.r")

geneLists <- c("Normal_skin_pos", "BCC","GWAS_psoriasis","IL17","IL12_23","TNF","IL36_MyD88","IFNg","MHC_classI","TLR", "IBD_mucosa")
for(geneL in geneLists) {
  #assign(paste(geneL, "_genes",sep=""), read.table(paste("/Users/so11/phd/psoriasis/selection_analyses/gene_lists/", geneL, ".txt", sep="")))
  assign(paste(geneL, "_genes", sep=""), read.table(paste("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/08_selection_analyses/pathway_dnds_gene_lists/", geneL, ".txt", sep="")))
}

for(geneL in geneLists) {
  #assign(paste("dnds_",geneL,sep=""), dndscv_pvalspois(mutation_list[, c(1:5)], get(paste(geneL, "_genes", sep=""))$V1,max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf,refdb="/Users/so11/phd/psoriasis/resources/RefCDS_human_GRCh38.p12.rda", cv=NULL,))
    assign(paste("dnds_",geneL,sep=""), dndscv_pvalspois(mutation_list[, c(1:5)], get(paste(geneL, "_genes", sep=""))$V1,max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf,refdb=refcds_38, cv=scores))
}

d38 = dndscv(all_mutations, refdb=refcds_38, cv=scores,max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf) 
                                     
df <- data.frame()
for(geneL in geneLists) {
  assign(paste("sel_cv_", geneL, sep=""), get(paste("dnds_",geneL,sep=""))$sel_cv)
  
  df_tmp <- data.frame(get(paste("dnds_",geneL,sep=""))$globaldnds)
  df_tmp$geneList <- geneL
  df <- rbind(df, df_tmp)
}
mis_and_trunc <- df[df$name %in% c("wmis", "wtru"),]

P_mis <- as.numeric()
P_trunc <- as.numeric()

for(geneL in geneLists) {
  P_mis <- c(P_mis, summary(get(paste("dnds_", geneL, sep=""))$poissmodel)$coefficients[c("wmis"),4] )
  P_trunc <- c(P_trunc, summary(get(paste("dnds_", geneL, sep=""))$poissmodel)$coefficients[c("wtru"),4] )
}
mis_and_trunc$p_pois <- NA
mis_and_trunc$p_pois[seq(from=1, to=nrow(mis_and_trunc), by=2)] <- P_mis
mis_and_trunc$p_pois[seq(from=2, to=nrow(mis_and_trunc), by=2)] <- P_trunc
mis_and_trunc$q_BHcorr <- p.adjust(mis_and_trunc$p_pois, method="BH")

write.table(mis_and_trunc, "/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/pathway_dnds_mis_and_trunc_w_covs.txt", row.names=F, quote=F, sep="\t")














