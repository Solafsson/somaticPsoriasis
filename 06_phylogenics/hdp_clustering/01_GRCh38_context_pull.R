### Script to pull sequence context, gene ID, and transcript strand 
### Using GRCh38, replacing Context_pull_build37.pl script for running NDP

.libPaths("/lustre/scratch119/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")

library(data.table)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(Biostrings)

args = commandArgs(trailingOnly=TRUE)
hdp_dir=args[1]
patient <- args[2]


### Get context pull input files
pull.dir <- paste(hdp_dir, patient, "/", sep="")
pull.files <- list.files(pull.dir, full.names = T, pattern = paste(patient,"context_pull.txt", sep="_"))

### Loop through pull files, generate and write outputs
for (i in 1:length(pull.files)) {
  this.file = pull.files[i]
  this.pos.list = read.table(this.file, h=T)
  colnames(this.pos.list ) <- c("chrom", "pos", "ref", "alt")
  context_list <- this.pos.list  %>%
    dplyr::mutate(start = pos - 10,
                  end = pos + 10) %>%
    dplyr::select(chrom, start,end)
  this_range <- GenomicRanges::makeGRangesFromDataFrame(context_list)
  this_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, this_range)
  out_seqs <- GenomicRanges::as.data.frame(this_seq)
  out.table <- bind_cols(this.pos.list, out_seqs)
  colnames(out.table) <- c("chrom", "pos", "ref", "alt", "CONTEXT")
  write.table(out.table, file = paste0(pull.dir, patient, "_mut_context_GRCh38.txt"), quote = F, row.names = F, sep = "\t")
}
