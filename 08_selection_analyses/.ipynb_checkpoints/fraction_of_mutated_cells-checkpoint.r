
# Input files
exome_and_tgs_mutations = "bld_exometgs_subsindels_combined.tsv"
exome_mutations = "bld_exome_caveman_and_pindel_calls.tsv"
target_genes_file = "normal_baits_genes.txt"
bladder_cancer_genes_file = "Bladder_cancer_genes_PMID_28988769_29056346.tsv"
targeted_coverage_file = "bld_targeted_picard_coverage.tsv"
exome_coverage_file = "bld_exome_picard_coverage.tsv"
tcgablca_file = "mc3.v0.2.8.PUBLIC.uniqueperdonor.BLCA.6cols.txt"
patient_file = "bladder_patient_info_2019-10-30.csv"
lcm_file = "2019-10-01_LCM_database.rds"
# Loading mutations (substitutions and indels from TGS and WES)
muts.all = read.table(exome_and_tgs_mutations, sep="\t", header=1, stringsAsFactors=F)
targetgenes = c(read.table(target_genes_file, header=0, sep="\t", stringsAsFactors=F)[,1], "CDKN2A.p14arf", "CDKN2A.p16INK4a") # Adding the two CDKN2A isoforms
muts.exome = read.table(exome_mutations, sep="\t", header=1, stringsAsFactors=F)
# Known bladder cancer genes
bladder_cancer_genes = read.table(bladder_cancer_genes_file, header=0, sep="\t", stringsAsFactors=F)[,1]
# Selecting urothelial samples
patientdata = read.table(patient_file, header=1, sep=",", stringsAsFactors=F)
lcmdata = readRDS(lcm_file)
urot_ids = as.vector(lcmdata$SupplierSampleName[lcmdata$Feature=="Urothelium" & as.vector(lcmdata$Donor) %in% patientdata$external_id[patientdata$patient_type=="transplant"]])
mutations = muts.all[which(muts.all$sampleID %in% urot_ids),]
mutations$sampleID = substr(mutations$sampleID,1,7) # Using patient identifiers instead of sample identifiers
mutations = unique(mutations) # Unique mutations per patient



# c. % of mutant urothelium (we will use targeted sequenced samples with median coverage >=50 from transplant donors >=50 years)

min_coverage = 50 # Only considering samples with this minimum coverage
min_age = 50 # Only donors with this minimum age
coverage_table1 = read.table(targeted_coverage_file, header=1, sep="\t", stringsAsFactors=F)
coverage_table2 = read.table(exome_coverage_file, header=1, sep="\t", stringsAsFactors=F)
enough_cov = unique(c(coverage_table1$SAMPLE_ID[which(coverage_table1$MEDIAN_TARGET_COVERAGE>=min_coverage)],
                      coverage_table2$SAMPLE_ID[which(coverage_table2$MEDIAN_TARGET_COVERAGE>=min_coverage)]))

urot_ids = intersect(enough_cov, as.vector(lcmdata$SupplierSampleName[lcmdata$Feature=="Urothelium" & as.vector(lcmdata$Donor) %in% patientdata$external_id[patientdata$patient_type=="transplant"]]))
urot_ids = urot_ids[substr(urot_ids,1,7) %in% patientdata$internal_id[patientdata$age>=min_age]]
num_samples = sum((lcmdata$SentForTargeted=="Y" | lcmdata$SentForExome=="Y") & lcmdata$SupplierSampleName %in% urot_ids)

mutations = dndscv(muts.all[which(muts.all$sampleID %in% urot_ids),1:5], gene_list=targetgenes, max_muts_per_gene_per_sample=Inf, max_coding_muts_per_sample=Inf, outp=1)$annotmuts # Annotating the mutations
mutations$str = paste(mutations$sampleID, mutations$chr, mutations$pos, mutations$mut, sep=":")
vafs = setNames(muts.all$vaf, paste(muts.all$sampleID, muts.all$chr, muts.all$pos, muts.all$mut, sep=":"))
mutations$vaf = vafs[mutations$str]

# Correction for X-chr genes
xchr_genes = unique(dndsout_ref$annotmuts$gene[dndsout_ref$annotmuts$gene %in% genes2plot & dndsout_ref$annotmuts$chr=="X"])
male_patients = patientdata$internal_id[patientdata$gender=="Male"]

# Initialising analysis per gene
fraction_mutant_cells = array(0,dim=c(2,length(genes2plot)))
colnames(fraction_mutant_cells) = genes2plot
rownames(fraction_mutant_cells) = c("highbd","lowbd")
ns = mutations[which(mutations$impact!="Synonymous" & mutations$gene %in% genes2plot),] # Table with all non-synonymous mutations

# Initialising analysis per sample per gene (for the calculation of totals across genes)
fraction_mutant_cells_persample_low = array(0, dim=c(num_samples,length(genes2plot)))
s = unique(mutations$sampleID)
rownames(fraction_mutant_cells_persample_low) = c(s, rep("",num_samples-length(s)))
colnames(fraction_mutant_cells_persample_low) = genes2plot
fraction_mutant_cells_persample_high = fraction_mutant_cells_persample_low

for (j in 1:length(genes2plot)) {
  nsj = ns[ns$gene==genes2plot[j],]
  nspersample = split(nsj, f=nsj$sampleID)
  if (nrow(nsj)==0) {
    fraction_mutant_cells[,j] = 0
  } else {
    if (!(genes2plot[j] %in% xchr_genes)) { # Genes in autosomes
      
      aux = sapply(nspersample, function(x) c(min(1,sum(x$vaf*2)), min(1,sum(x$vaf))))
      fraction_mutant_cells[,j] = rowSums(aux) / num_samples
      fraction_mutant_cells_persample_low[colnames(aux), genes2plot[j]] = aux[2,]
      fraction_mutant_cells_persample_high[colnames(aux), genes2plot[j]] = aux[1,]
      
    } else { # Genes in X-chr (different calculation for male and female patients)
      f = array(0,dim=c(length(nspersample),2)) # Initialise
      rownames(f) = names(nspersample)
      for (p in 1:length(nspersample)) {
        if (substr(nspersample[[p]]$sampleID[1],1,7) %in% male_patients) { # Male
          f[p,] = min(1,sum(nspersample[[p]]$vaf))
        } else { # Female (standard diploid case)
          f[p,] = c(min(1,sum(nspersample[[p]]$vaf*2)), min(1,sum(nspersample[[p]]$vaf)))
        }
      }
      fraction_mutant_cells[,j] = colSums(f) / num_samples
      fraction_mutant_cells_persample_low[rownames(f), genes2plot[j]] = f[,2]
      fraction_mutant_cells_persample_high[rownames(f), genes2plot[j]] = f[,1]
    }
  }
}
aux = rbind(fraction_mutant_cells[2,],fraction_mutant_cells[1,]-fraction_mutant_cells[2,]) * 100
pos = barplot(aux, las=2, col=c("white","indianred3"), border=NA, ylab="% mutant epithelium", ylim=c(0,6))
abline(v=verticalbar)