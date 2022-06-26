.libPaths("/lustre/scratch119/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")
library(MutationalPatterns)
library(BSgenome)
#Import Reference Genome
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(stringr)

## Define variables
binom_filter_dir <- "/lustre/scratch119/humgen/projects/psoriasis/binomial_filters/"
patients <- read.table("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/hdp/snv/patient_list.txt")
output_dir <- "/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/hdp/indels/"



indels <- data.frame()
for(pat in patients$V1) {
  pat_indels <- read.table(paste(binom_filter_dir, pat, "/", pat, "_genotype_indels.txt", sep=""), h=T)
  pat_indels$patientID <- pat
  indels <- rbind(indels, data.frame(rownames(pat_indels), pat_indels$patientID))
}

colnames(indels) <- c("mutationID", "patientID")

indels$Chr <- unlist(strsplit(indels$mutationID, split=":"))[c(T,F,F,F)]
indels$Pos <- unlist(strsplit(indels$mutationID, split=":"))[c(F,T,F,F)]
indels$Ref <- unlist(strsplit(indels$mutationID, split=":"))[c(F,F,T,F)]
indels$Alt <- unlist(strsplit(indels$mutationID, split=":"))[c(F,F,F,T)]


Grange_branches <- makeGRangesListFromDataFrame(indels, split.field = "patientID", keep.extra.columns = T, ignore.strand = T, seqnames.field = "Chr", 
                                                start.field = "Pos", end.field = "Pos")
GenomeInfoDb::genome(Grange_branches) = 'hg38'

indel_grl <- get_indel_context(Grange_branches, ref_genome)
indel_counts <- count_indel_contexts(indel_grl)
indel_counts <- data.frame(indel_counts)

## Combine the count for PUVA vs non-PUVA exposed patients
clone_sigs <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/clone_mutation_burden.txt", h=T)
x <- clone_sigs[clone_sigs$Psoralens>100,]
puva_exp_patients <- unique(unlist(strsplit(x$CloneID, spli="_"))[c(T,F)])

exp <- indel_counts[,colnames(indel_counts) %in% puva_exp_patients]
nexp <- indel_counts[,!(colnames(indel_counts) %in% puva_exp_patients)]

sums <- data.frame(Exposed=rowSums(exp), NonExposed=rowSums(nexp))
plot_indel_contexts(sums, condensed = TRUE)
