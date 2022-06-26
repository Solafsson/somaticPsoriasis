
## The purpose of this script is to estimate the sensitivity of the mutation calls
## From the technical duplicates
## Assuming the same sensitivity in both samples, the maximum likelihood estimate for 
## the sensitivity is:
## S = 2*n2/(n1 + 2*n2)
## Where n2 is the number of mutations called in both samples and n1 is the sum 
## of mutations called in only one sample.
.libPaths("/lustre/scratch119/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")
library(cowplot)
library(ggplot2)

sample_meta <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/sample_info/sample_meta.txt", h=T)
dups <- sample_meta[sample_meta$excl_criteria=="Technical_duplicate" & !is.na(sample_meta$excl_criteria),]
dups$original <- gsub("_repl", "", dups$sampleID)
dups$original[dups$original=="P03L_T_10"] <- "P03L_T_5"
dups$original[dups$original=="P03L_T_9"] <- "P03L_T_4"

## P03L_T_10 is a duplicate of P03L_T_5 and P03L_T_9 is a duplicate of P03L_T_4
## They don't have the _repl names because I cut the same histological features 
## twice by mistake. 

cov <- read.table("/lustre/scratch119/humgen/projects/psoriasis/coverage_calculation/samtools_depth_coverage.txt")
cov <- cov[,c(1,7)]
colnames(cov) <- c("sampleID", "Median_Coverage")

S <- as.numeric()
cov_original <- as.numeric()
cov_dup <- as.numeric()
lower <- as.numeric()
for(i in 1:nrow(dups)) {
  patient <- dups$patient_ID[i]
  VAF_matrix <- read.table(paste("/lustre/scratch119/humgen/projects/psoriasis/binomial_filters/", patient, "/",patient,"_VAFs_pass_allMuts.txt", sep=""), h=T)
  genotype_bin <- read.table(paste("/lustre/scratch119/humgen/projects/psoriasis/binomial_filters/", patient, "/",patient,"_genotype_allMuts.txt", sep=""), h=T)
  
  theseTwo_vaf <- VAF_matrix[,c(dups$sampleID[i], dups$original[i])]
  colnames(theseTwo_vaf) <- c("Duplicate", "Original")
  theseTwo_bin <- genotype_bin[,c(dups$sampleID[i], dups$original[i])]
  colnames(theseTwo_bin) <- c("Duplicate", "Original")
  
  x <- theseTwo_bin[!(theseTwo_bin$Duplicate==0 & theseTwo_bin$Original==0),]
  n2 <- nrow(x[x$Duplicate>0 & x$Original>0,])
  n1 <- nrow(x[(x$Duplicate==0 & x$Original==1) | (x$Duplicate==1 & x$Original==0),])
  S[i] = 2*n2/(n1 + 2*n2)
  cov_original[i] <- cov$Median_Coverage[cov$sampleID==dups$original[i]]
  cov_dup[i] <- cov$Median_Coverage[cov$sampleID==dups$sampleID[i]]
  if(cov_original[i]<cov_dup[i]) {
    lower[i] <- cov_original[i]
  } else {
    lower[i] <- cov_dup[i]
  }
}

## 1 cut to the original, next to the duplicate, then to the original etc.
## as opposed to first cut all samples for original and then all for duplicate
braid_sampling <- c(0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0)
braid_sampling <- ifelse(braid_sampling==0, "No", "Yes")

sensitivity <- data.frame(Original=dups$original, cov_original, Duplicate=dups$sampleID, cov_dup,lower, S, braid_sampling)
write.table(sensitivity, file="/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/Sensitivity_technical_duplicates.txt", sep="\t", row.names=F)

ggplot(sensitivity, aes(x=lower,y=S, colour=braid_sampling )) + geom_point() + ylim(0,1) + 
  theme_bw(base_size = 14) + labs(x="Lower coverage of pair", y="Sensitivity", colour="Braid sampling:") + 
  theme(legend.position = "top") + geom_vline(xintercept=median(cov$Median_Coverage)) 


ggplot(theseTwo, aes(x=Original, y=Duplicate)) + geom_point(alpha=0.5) + theme_classic() + 
  labs(x="Original VAF", y="Technical duplicate VAF")













