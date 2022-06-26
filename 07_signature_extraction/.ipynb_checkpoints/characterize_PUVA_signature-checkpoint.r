
## The purpose of this script is to characterize the mutational signature
## associated with PUVA treatment. I will use the whole-genome sequencing
## data from patients 18, 21 and 34 for this. 

## Relying heavily on this tutorial of the MutationalPatterns package:
## https://bioconductor.org/packages/release/bioc/vignettes/MutationalPatterns/inst/doc/Introduction_to_MutationalPatterns.html#transcriptional-strand-bias-analysis


.libPaths("/lustre/scratch119/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")

library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(grid)
library(gridExtra)
library(MutationalPatterns)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
options(stringsAsFactors = F)


## Define a few variables
genes_hg38 <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
sample_meta <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/sample_info/sample_meta_wgs.txt", h=T, stringsAsFactors = F)
binomial_dir="/lustre/scratch119/humgen/projects/psoriasis/binomial_filters/"
repl_timing_dir="/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/07_signature_extraction/replication_timing/"
puva_dir="/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/"

sample_names <- sample_meta$sampleID
patient_list <- sample_meta$patient_ID

## Start by converting the binary genotype matrices into vcfs
# patient <- "patient18_WGS"

for(patient in unique(patient_list)) {
  
  genotype_bin <- read.table(paste(binomial_dir, patient, "/", patient, "_genotype_allMuts.txt", sep=""), h=T)
  
  for(sample in colnames(genotype_bin)) {
    
    sample_muts <- rownames(genotype_bin[genotype_bin[,which(colnames(genotype_bin)==sample)]==1,])
    sample_muts <- sample_muts[!is.na(sample_muts)]
    
    sample_df <- data.frame(chr=unlist(strsplit(sample_muts, split=":"))[c(T,F,F,F)],
                            start=unlist(strsplit(sample_muts, split=":"))[c(F,T,F,F)],
                            end=unlist(strsplit(sample_muts, split=":"))[c(F,T,F,F)],
                            REF=unlist(strsplit(sample_muts, split=":"))[c(F,F,T,F)],
                            ALT=unlist(strsplit(sample_muts, split=":"))[c(F,F,F,T)])
    
    sample_df$QUAL <- "."
    sample_df$FILTER <- "PASS"
    sample_df$INFO <- "."
    sample_df$FORMAT <- "."
    sample_df$NORMAL <- "."
    sample_df$TUMOUR <- "."
    colnames(sample_df)[colnames(sample_df)=="end"] <- "ID"
    sample_df <- sample_df[nchar(sample_df$REF)==1,]
    
    write.table(sample_df, file=paste(puva_dir, sample, "_vcf_no_header.txt", sep=""), row.names = F, col.names = F, quote = F, sep="\t")
    system(paste("cat ", puva_dir, "vcf_header.txt ", puva_dir, sample, "_vcf_no_header.txt"," > ", puva_dir,
                 sample, ".vcf", sep=""))
  }
}

vcf_files <- paste(puva_dir, sample_names, ".vcf", sep="")

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome, predefined_dbs_mbs=T)


## Plot the mutational profiles of PUVA exposed samples:
###########################################################

type_occurrences <- mut_type_occurrences(grl, ref_genome)
plot_spectrum(type_occurrences, CT = TRUE, 
              indv_points = TRUE, legend = T, by=patient_list)

mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
plot_96_profile(mut_mat[,c(1:4)])
plot_96_profile(mut_mat[,c(5:8)])
plot_96_profile(mut_mat[,c(9:16)])

## Look at larger context than the trinucleotide context
mut_mat_ext_context <- mut_matrix(grl, ref_genome, extension = 2)
plot_profile_heatmap(mut_mat_ext_context, by=patient_list)

## We could in principle look beyond the pentanucleotide model but
## this quickly gets very messy... 
mut_mat_ext_context3 <- mut_matrix(grl, ref_genome, extension = 3)
plot_profile_heatmap(mut_mat_ext_context3)

## A second way of showing the effect of sequence context is to plot a 'river plot.
plot_river(mut_mat_ext_context)
plot_river(mut_mat_ext_context[,c(1:4)])
plot_river(mut_mat_ext_context[,c(5:8)])
plot_river(mut_mat_ext_context[,c(9:16)])


## Indels
#####################
## We can see that the indel mutation spectrum is dominated by 
## COSMIC indel signature 1, attributed to UV-damage. There is 
## no effect of PUVA treatment on the indel signature. 
indel_grl <- read_vcfs_as_granges(vcf_files, sample_names, 
                                  ref_genome, type = "indel")
indel_grl <- get_indel_context(indel_grl, ref_genome)
indel_counts <- count_indel_contexts(indel_grl)

plot_indel_contexts(indel_counts[,c(1:4)], condensed = F)+ theme(legend.position="top")
plot_indel_contexts(indel_counts[,c(5:8)], condensed = TRUE)
plot_indel_contexts(indel_counts[,c(9:16)], condensed = TRUE)
plot_main_indel_contexts(indel_counts)

## DBS mutations
#####################
dbs_grl <- get_mut_type(grl, type = "dbs")

predefined_dbs_grl <- read_vcfs_as_granges(vcf_files, sample_names, 
                                           ref_genome, type = "dbs",
                                           predefined_dbs_mbs = TRUE)

dbs_grl <- get_dbs_context(predefined_dbs_grl)



## Use MutationalPatterns to carry out a basic NMF-based
## signature extraction

signatures = get_known_signatures(genome = "GRCh38")

library("NMF")
mut_mat <- mut_mat + 0.0001
estimate <- nmf(mut_mat, rank = 2:5, method = "brunet", 
                nrun = 10, seed = 123456, .opt = "v-p")

plot(estimate)

nmf_res <- extract_signatures(mut_mat, rank = 4, nrun = 10, single_core = TRUE)
colnames(nmf_res$signatures) <- c("Signature A", "Signature B")
rownames(nmf_res$contribution) <- c("Signature A", "Signature B")
nmf_res <- rename_nmf_signatures(nmf_res, signatures, cutoff = 0.85)
colnames(nmf_res$signatures)
plot_96_profile(nmf_res$signatures, condensed = TRUE)


## The above is absolute trash. 
## Try fitting the mutation matrix to the COSMIC mutational signatures:
hdp_puva_component <- read.table("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/PUVA_signature_hdp_component.txt")
colnames(hdp_puva_component) <- "PUVA"
signatures <- as.matrix(cbind(signatures, hdp_puva_component))


fit_res <- fit_to_signatures(mut_mat, signatures[,c("SBS1", "SBS5", "SBS7a", "SBS7b",
                                                    "SBS7c","SBS7d", "PUVA")])

plot_contribution(fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
) + theme(axis.text.x = element_text(angle=90))

plot_contribution(fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "relative"
) + theme(axis.text.x = element_text(angle=90))

plot_original_vs_reconstructed(mut_mat, fit_res$reconstructed, 
                               y_intercept = 0.95)

strict_refit <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.004)

fit_res_strict <- strict_refit$fit_res
plot_contribution(fit_res_strict$contribution,
                  coord_flip = FALSE,
                  mode = "absolute"
)

## SBS93 is interesting. Try including that as a prior in the 
## HDP extraction process. Found in ESCCs. Could explain that
## one HDP component that was strange in the hdp extraction.
## see https://www.biorxiv.org/content/10.1101/2020.12.13.422570v1.full.pdf


## Signature potential damage analysis
##############################################
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## Entrez gene IDs for c(TP53, NOTCH1,NOTCH2, PPM1D,BRAF)
gene_ids <- c(7157, 4851,4853,8493,673)

contexts <- rownames(mut_mat)
context_mismatches <- context_potential_damage_analysis(contexts, txdb, ref_genome, gene_ids)
head(context_mismatches)
sig_damage <- signature_potential_damage_analysis(signatures[, c("SBS7a", "SBS7b", "PUVA")], contexts, context_mismatches)
# A normalized ratio of 2 for "stop gain" mutations, means that a signature 
# is twice as likely to cause "stop gain" mutations, compared to a completely 
# random "flat" signature. 
sig_damage[1:4,c(1,4)]


## Strand bias analyses
##########################


## First look for a transcriptional strand bias
genes_hg38 <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
mut_mat_s <- mut_matrix_stranded(grl, ref_genome, genes_hg38)
plot_192_profile(mut_mat_s[,1:4]) + theme_classic(base_size = 14) + theme(legend.position="top", axis.text.x = element_text(angle=90, size=8))
plot_192_profile(mut_mat_s[,5:8]) + theme_classic(base_size = 14) + theme(legend.position="top", axis.text.x = element_text(angle=90, size=8))
plot_192_profile(mut_mat_s[,16]) + theme_classic(base_size = 14) + theme(legend.position="top", axis.text.x = element_text(angle=90, size=8))



strand_counts <- strand_occurrences(mut_mat_s)
strand_bias <- strand_bias_test(strand_counts)
strand_bias

plot_strand(strand_counts, mode = "relative")
plot_strand_bias(strand_bias, sig_type = "p")

## Next look for replication strand bias. 
## Use the file from Fede. He got it from here:
## http://www.cell.com/abstract/S0092-8674(15)01714-6

## In this context, left is synonymous to the leading strand
## because we are always using the reference as a point of reference.
## I observe an enrichment of PUVA related mutations on the leading (left)
## strand. 

repli_timing <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/07_signature_extraction/strand_asymmetries/replication_timing.repdirhg38.txt")
repli_timing$Class <- NA
repli_timing$Class[repli_timing$V5==1] <- "left"
repli_timing$Class[repli_timing$V6==1] <- "right"

repli_strand <- repli_timing[!is.na(repli_timing$Class), c("V1", "V2", "V3", "Class")]
repli_strand$Ratio <- 0
colnames(repli_strand)=c("Chr", "Start", "Stop", "Class", "Ratio")

# Store in GRanges object
repli_strand_granges <- GRanges(
  seqnames = repli_strand$Chr,
  ranges = IRanges(
    start = repli_strand$Start + 1,
    end = repli_strand$Stop
  ),
  strand_info = repli_strand$Class
)
# UCSC seqlevelsstyle
seqlevelsStyle(repli_strand_granges) <- "UCSC"
repli_strand_granges

repli_strand_granges$strand_info <- factor(repli_strand_granges$strand_info,
                                           levels = c("right", "left")
)


mut_mat_s_rep <- mut_matrix_stranded(grl, ref_genome, repli_strand_granges,
                                     mode = "replication"
)
strand_counts_rep <- strand_occurrences(mut_mat_s_rep)
strand_bias_rep <- strand_bias_test(strand_counts_rep)

ps1 <- plot_strand(strand_counts_rep, mode = "relative")
ps2 <- plot_strand_bias(strand_bias_rep)
grid.arrange(ps1, ps2)

## Test this out with the replication origin file that comes
## with the MutationalPatterns package. It is presumably in hg19
## so needs lifing over... 
repli_file <- system.file("extdata/ReplicationDirectionRegions.bed",
                          package = "MutationalPatterns"
)
repli_strand <- read.table(repli_file, header = TRUE)

write.table(repli_strand, 
            file="/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/07_signature_extraction/strand_asymmetries/MutationalPatterns_repli_strand.txt", sep="\t", row.names = F, quote = F)

# awk 'NR>1 {OFS="\t"; print "chr"$1, $2-1, $3-1, $4, $5}' < MutationalPatterns_repli_strand.txt > MutationalPatterns_repli_strandhg37.bed
# /software/team152/liftover/liftOver MutationalPatterns_repli_strandhg37.bed /software/team152/liftover/hg19ToHg38.over.chain.gz MutationalPatterns_repli_strandhg38_lifted.bed MutationalPatterns_repli_strandhg38_unlifted.bed
# awk 'BEGIN {OFS="\t"; print "Chr", "Start", "Stop", "Class", "Ratio"} {OFS="\t"; print $1, $2+1, $3+1, $4, $5}' < MutationalPatterns_repli_strandhg38_lifted.bed > MutationalPatterns_repli_strandhg38.txt

repli_strand <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/07_signature_extraction/strand_asymmetries/MutationalPatterns_repli_strandhg38.txt", h=T)

# Store in GRanges object
repli_strand_granges <- GRanges(
  seqnames = repli_strand$Chr,
  ranges = IRanges(
    start = repli_strand$Start + 1,
    end = repli_strand$Stop
  ),
  strand_info = repli_strand$Class
)
# UCSC seqlevelsstyle
seqlevelsStyle(repli_strand_granges) <- "UCSC"
repli_strand_granges

repli_strand_granges$strand_info <- factor(repli_strand_granges$strand_info,
                                           levels = c("right", "left")
)


mut_mat_s_rep <- mut_matrix_stranded(grl, ref_genome, repli_strand_granges,
                                     mode = "replication"
)
strand_counts_rep <- strand_occurrences(mut_mat_s_rep)
strand_bias_rep <- strand_bias_test(strand_counts_rep)
strand_bias_rep

ps1 <- plot_strand(strand_counts_rep, mode = "relative")
ps2 <- plot_strand_bias(strand_bias_rep)
grid.arrange(ps1, ps2)


## Replication timing analysis
####################################

## Read in the results from puva_replication_timing.sh

rt_results <- data.frame()
for(sample in sample_names) {
  rt <- read.table(paste(puva_dir, sample, "_RT_results.txt", sep=""))
  rt$sample <- sample
  rt_results <- rbind(rt_results, rt)
}

ggplot(rt_results, aes(x=V1, y=log10(V2))) + geom_boxplot() +geom_point() +
  geom_line(aes(group=sample, colour=sample)) + theme_classic(base_size = 14) +
  labs(y="log10(Mutations per TA/AT site)", x="")


## Save selected objects for plotting in make_manuscript_figures.r
save(grl, mut_mat_ext_context,mut_mat_s, strand_counts, mut_mat_s_rep,strand_counts_rep, rt_results, 
     file = "/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/puva_characterization.data.RData")


## Look for evidence of transcription coupled damage
######################################################
## See figure 7A from http://www.cell.com/abstract/S0092-8674(15)01714-6

# Compile a list of genes highly expressed in the skin. Include strand information.
# Create 1kb bins 10kb up and downstream of the TSSs
# Count the number of T>A, T>C and T>G mutations at TpA sites on the transcribed strand
# and the number of A>T, A>G and A>C mutations at ApT sites on the untranscribed strand
# if the gene is on the (+) strand, the transcribed strand is the complement of the 
# reference strand. 

## First read in GTEx results and merge with gencode
gtex <- read.csv("/lustre/scratch119/humgen/projects/psoriasis/resources/GTEx_skin_expression.csv", h=T)
colnames(gtex) <- c("TranscriptID", "Gene", "Skin_no_sun", "Skin_sun")

tss<-read.table("/lustre/scratch119/humgen/projects/psoriasis/resources/gencode.v27.annotation_chr_pos_strand_geneid_gene_name.gtf",h=F)
colnames(tss)<-c("Chromosome_phe","START_GENE","END_GENE","STRAND","Phenotype_ID","HGNC")
tss$TSS<-ifelse(tss$STRAND=="+",tss$START_GENE-1,tss$END_GENE-1)

## Extract the location and TSS for the 10% of genes with highest expression
## in the skin
comb <- merge(gtex, tss, by.x="Gene", by.y="HGNC")
comb <- comb[order(comb$Skin_sun, decreasing = T),]
comb <- comb[comb$Chromosome_phe!="chrM",]
d <- duplicated(comb$Gene)
comb <- comb[!d,]  ## Keeps the highest expressed transcript

# Use the genes in the top 10% of expressed transcripts in GTEx
toUse <- comb[1:(floor(nrow(comb)/10)),]

## For each of these genes, I want to write 20 entries to a bed file,
## ten 1kb segments upstream and downstream of the TSS. 
chr_vector <- as.character()
start <- as.numeric()
end <- as.numeric()
tss <- as.character()
strand <- as.character()
label <- as.character()
t <- 1
for(i in 1:nrow(toUse)) {
  

  if(toUse$STRAND[i]=="-") {
    # Create 1kb bins upstream of the TSS
    for(j in 1:10) {
      chr_vector[t] <- toUse$Chromosome_phe[i]
      start[t] <- toUse$END_GENE[i] + (j-1)*1000
      end[t] <- toUse$END_GENE[i] + (j)*1000
      tss[t] <- "upstream"
      strand[t] <- "-"
      label[t] <- paste("-",j,"kb", sep="")
      t <- t+1
    }
    # Create 1kb bins downstream of the TSS
    for(k in 1:10) {
      chr_vector[t] <- toUse$Chromosome_phe[i]
      start[t] <- toUse$END_GENE[i] - (k)*1000
      end[t] <- toUse$END_GENE[i] - (k-1)*1000
      tss[t] <- "downstream"
      strand[t] <- "-"
      label[t] <- paste("+",k,"kb", sep="")
      t <- t+1
    }
  } else if(toUse$STRAND[i]=="+") {
    # Create 1kb bins upstream of the TSS
    for(j in 1:10) {
      chr_vector[t] <- toUse$Chromosome_phe[i]
      start[t] <- toUse$START_GENE[i] - (j)*1000
      end[t] <- toUse$START_GENE[i] - (j-1)*1000
      tss[t] <- "upstream"
      strand[t] <- "+"
      label[t] <- paste("-",j,"kb", sep="")
      t <- t+1
    }
    # Create 1kb bins downstream of the TSS
    for(k in 1:10) {
      chr_vector[t] <- toUse$Chromosome_phe[i]
      start[t] <- toUse$START_GENE[i] + (k-1)*1000
      end[t] <- toUse$START_GENE[i] + (k)*1000
      tss[t] <- "downstream"
      strand[t] <- "+"
      label[t] <- paste("+",k,"kb", sep="")
      t <- t+1
    }
  } else {
    stop("Strand info missing")
  }
}

bed_file <- data.frame(chr_vector, start, end, tss, strand, label)

write.table(bed_file, 
            file = "/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/top10pc_expr_genes_skin.bed", 
            sep="\t", quote = F, row.names = F, col.names = F)

