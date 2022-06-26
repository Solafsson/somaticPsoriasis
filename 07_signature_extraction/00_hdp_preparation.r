## This script reads in mutations assigned to branches of a phylogenic tree in /nfs/users/nfs_s/so11/phd/somatic_ibd_p1/phylogenics/assign_mutations_to_branches.r
## It creates a mutational matrix on the correct form for subsequent HDP signature extraction. 
## As hdp_posterior_sampling_worker.r is run in a job array, this script prevents each job from having to 
## repeat the generation of the mutational matrix. Data can instead be loaded, which should be faster. 

## I only include branches that have more than <branch_threshold> mutations. 


.libPaths("/lustre/scratch119/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")
library(MutationalPatterns)
library(BSgenome)
#Import Reference Genome
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(stringr)
# see package vignettes with
# browseVignettes('hdp')

## Define variables
branch_threshold <- 20
hdp_directory <- "/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/"


patients <- read.table("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/hdp/snv/patient_list.txt")
output_dir <- "/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/hdp/snv/"


branches <- data.frame()
for(pat in patients$V1) {
  branch <- read.table(paste(hdp_directory, pat, "/ndp_", pat, "_2021_12_09sbs_assignment.txt", sep=""), h=T)
  branches <- rbind(branches, branch)
}

use_branches <- as.character(subset(data.frame(table(branches$ClusterID)), Freq > branch_threshold)[,1])
branches2 <- subset(branches, as.character(ClusterID) %in% use_branches & !is.na(Cluster))


#branch_data <- as.data.frame(str_split_fixed(branches2$mutationID, "_", n=4))
branch_data <- branches2[, c("Chr", "Pos", "Ref", "Alt")]
colnames(branch_data) <- c("Chr", "Pos", "REF", "ALT")
branch_data$Pos <- as.numeric(as.character(branch_data$Pos))
branch_data$Cluster <- branches2$ClusterID
#branch_data$Chr <- paste("chr", branch_data$Chr, sep="")

Grange_branches <- makeGRangesListFromDataFrame(branch_data, split.field = "Cluster", keep.extra.columns = T, ignore.strand = T, seqnames.field = "Chr", 
                                                start.field = "Pos", end.field = "Pos")
GenomeInfoDb::genome(Grange_branches) = 'hg38'

#sub_types <- mut_type_occurrences(Grange_branches, ref_genome)
mut_mat <- mut_matrix(Grange_branches, ref_genome)
mut_mat <- mut_mat[, colSums(mut_mat)!=0]
mut_mat_t <- t(mut_mat)



save(mut_mat_t, file=paste(output_dir, "branch_mutational_matrix.txt", sep=""))

