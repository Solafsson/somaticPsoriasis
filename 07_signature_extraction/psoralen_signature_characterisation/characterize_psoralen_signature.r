
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
plot_96_profile(mut_mat[,c(5:8)], condensed = T)
plot_96_profile(mut_mat[,c(9:16)],condensed = T)

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

dbs_grl <- get_dbs_context(dbs_grl)
dbs_counts <- count_dbs_contexts(dbs_grl)
plot_dbs_contexts(dbs_counts, same_y = TRUE)


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

signatures = get_known_signatures(genome = "GRCh38")
hdp_puva_component <- read.table("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/PUVA_signature_hdp_component.txt")
colnames(hdp_puva_component) <- "PUVA"
signatures <- as.matrix(cbind(signatures, hdp_puva_component))

## Entrez gene IDs for c(TP53, NOTCH1,NOTCH2, PPM1D,BRAF)
gene_ids <- c(7157, 4851,4853,8493,673)

## COSMIC census genes:
gene_ids <- c(29974,10006,25,27,57007,2181,23305,90,92,4301,4299,3899,27125,10142,207,208,10000,217,238,139285,286,324,9582,367,369,23092,2909,394,9639,55160,23365,8289,57492,196528,405,79058,171023,55252,466,471,472,476,492,545,546,8312,8313,567,8314,580,581,11177,8915,53335,64919,596,83596,602,604,605,607,283149,9774,54880,63035,613,330,57448,641,653,657,673,672,675,8019,23476,83990,694,7832,695,701,145788,776,811,23261,124583,84433,833,836,841,842,863,865,867,868,23624,8030,57820,892,595,894,896,898,1233,1236,30835,29126,940,972,973,974,79577,999,1008,1009,1015,51755,1019,1021,1026,1027,1029,1031,1045,1050,84902,79145,1106,1108,11200,26511,50515,23152,4261,6249,10978,1213,8218,168975,7555,4849,26047,11064,1277,1280,1281,1345,22849,1385,90993,64764,1387,64109,51340,23373,64784,1436,1441,114788,10664,1496,1499,1500,1501,8452,1523,7852,1540,1558,57105,1616,340578,1630,1639,1643,1649,4921,1662,1654,1655,1656,7913,54487,23405,3337,1785,1788,29102,1879,345930,8726,1956,107984923,3646,1974,1999,2000,2005,8178,2006,27436,2033,2034,2042,2045,2060,2064,2065,2066,23085,2068,2071,2072,2073,2078,2099,55500,2115,2118,2119,2120,2130,2131,2132,2146,7430,9715,51059,442444,2175,2176,2177,2178,2188,2189,355,2195,120114,79633,2199,80204,55294,2213,83417,2237,2242,54738,2260,11116,2263,2261,2264,2271,2272,81608,11328,201163,2313,2316,2322,2324,23048,3169,668,2308,2309,4303,27086,283150,10272,8880,2521,8522,2623,2624,2625,2735,8833,2767,2776,2778,9950,57120,2719,2262,10243,2903,2913,3020,3021,9709,23462,3091,3092,8358,8294,3105,3131,3159,8091,6927,3181,84376,3207,3209,3205,3227,3229,3237,3239,3265,3320,3326,3399,3417,3418,10644,3492,50802,3535,3551,10320,22806,3558,50615,3572,3575,3662,8471,91464,3685,3702,3716,3717,3718,221895,3725,7994,23522,11143,3762,5927,8242,7403,3791,2531,9817,57670,3799,3815,9314,1316,3817,4297,58508,8085,57082,90417,3845,3895,23185,3927,9113,26524,3932,3936,51176,23484,10186,3977,4000,4004,4005,4026,121227,53353,26065,4066,4067,8216,346389,4094,9935,10892,84441,5604,5605,6416,4214,9175,5594,4149,151963,4193,4194,2122,9968,4221,4233,4255,4286,4291,4292,4298,8028,10962,4300,4302,4330,3110,4352,57591,4436,2956,124540,4478,4515,2475,4582,94025,4585,4595,4602,4609,4610,4613,4615,4629,4627,4644,4654,55728,4665,4666,26960,4683,51517,8648,10499,8031,9611,9612,10397,4763,4771,4773,4780,4781,4791,4794,51199,7080,4841,4851,4853,4869,8013,4893,3084,64324,7468,54904,22978,4913,4914,4916,4926,8021,4928,256646,729262,728130,10215,4958,286530,26986,5049,79728,23598,5077,5079,5081,7849,55193,5087,5093,5108,80380,9659,5155,5156,5159,5187,84295,8929,8301,5290,5291,5295,5292,5324,5335,5371,5378,5395,5424,5426,5428,10721,25913,5450,5460,5468,8496,8493,5518,5537,5546,639,63976,7799,80243,5551,5566,5573,5579,25766,5396,11168,5727,5728,5753,5781,5783,5777,5787,5788,5789,5796,11122,114825,9444,9135,5879,5884,5885,5890,5894,5900,5903,5910,5914,5925,8241,64783,9401,5966,5979,55159,653489,6000,387,399,116028,57674,54894,6092,6098,6134,6146,6125,6184,340419,84870,861,862,6278,57167,51119,6385,6389,54949,6390,6391,6392,5413,23157,10801,6418,26040,23067,29072,9869,23451,6421,6424,6446,10019,6455,57698,140885,6495,10736,6497,10568,85414,4087,4088,4089,6597,6598,6602,6605,8243,6608,27044,92017,8651,6657,11166,92521,23013,8405,6714,9901,6427,6428,6760,26039,6756,6757,6759,10274,10735,6774,6777,6778,6491,6794,6801,51684,23512,6850,8148,6886,6887,79718,6926,6917,6938,6929,6934,8115,7006,54855,7015,80312,54790,7030,7942,10342,29844,7037,7048,9967,3195,30012,55654,7113,7114,3371,7128,8764,608,7150,7157,8626,7170,7171,7175,6955,84231,6957,6964,8805,5987,51592,9321,8295,7248,7249,7253,7307,51366,84101,9098,9101,8239,7409,7428,143187,7454,80304,11197,65268,7486,7490,25937,7507,7508,7514,7531,7704,55596,6935,463,7750,9203,55422,171017,353088,90827,25925,84133,8233)

## Genes mutated in SCCs and Normal skin:
gene_ids <- c(196528,171023,841,1029,2195,8085,4851,4853,8493,5925,8241,7157,8626,84962,4854)

contexts <- rownames(mut_mat)
context_mismatches <- context_potential_damage_analysis(contexts, txdb, ref_genome, gene_ids)
head(context_mismatches)
sig_damage <- signature_potential_damage_analysis(signatures[, c("SBS7b","PUVA")], contexts, context_mismatches)
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

plot_strand(strand_counts, mode = "absolute")
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


### READ IN EXPRESSION DATA FROM THE SKIN
##############################################

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


## Check PUVA mutation rate by expression level
######################################################

## Create expression bins to check the mutation rate in bins of ascending
## expression
bins <- comb[,c("Chromosome_phe", "START_GENE", "END_GENE", "Gene", "Skin_sun")]

sel_cv <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/sel_cv_dNdS_results.txt", h=T)
bins <- bins[bins$Gene %in% sel_cv$gene_name,]  ## Only consider protein coding genes. 

bins$nr <- 1:nrow(bins)
bins$expr_bin <- cut(bins$nr, breaks=10)
bins$expr_bin <- factor(bins$expr_bin, labels=c("Bin10","Bin9","Bin8","Bin7","Bin6","Bin5","Bin4","Bin3","Bin2","Bin1"))
bins$nr <- NULL
write.table(bins, file="/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/GTex_expression_bins.bed", sep="\t", col.names = F, row.names = F, quote = F)

## Run PUVA_expression.sh and read in the results
expression_res <- read.table("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/Mutation_rate_pr_expr_bin.txt", h=T)
expression_res$scaled <- expression_res$Mutation_Rate/expression_res$Mutation_Rate[1]

ggplot(expression_res, aes(x=Bin, y=scaled)) + geom_bar(position="dodge",stat="identity", fill="#E24E1B", colour="black") + 
  scale_y_continuous(expand=c(0,0), limits = c(0,1.2)) + theme_bw(base_size = 14) + 
  labs(y="Scaled mutation rate \n at TpA sites", x="Genes binned by increasing expression") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_hline(yintercept = 1, linetype=2)

## Look for evidence of transcription coupled damage
######################################################
## See figure 7A from http://www.cell.com/abstract/S0092-8674(15)01714-6

# Compile a list of genes highly expressed in the skin. Include strand information.
# Create 1kb bins 10kb up and downstream of the TSSs
# Count the number of T>A, T>C and T>G mutations at TpA sites on the transcribed strand
# and the number of A>T, A>G and A>C mutations at ApT sites on the untranscribed strand
# if the gene is on the (+) strand, the transcribed strand is the complement of the 
# reference strand. 
expr_quintile1 <-comb[comb$Gene %in% bins$Gene[bins$expr_bin=="Bin1" | bins$expr_bin=="Bin2"],]
expr_quintile2 <-comb[comb$Gene %in% bins$Gene[bins$expr_bin=="Bin3" | bins$expr_bin=="Bin4"],]
expr_quintile3 <-comb[comb$Gene %in% bins$Gene[bins$expr_bin=="Bin5" | bins$expr_bin=="Bin6"],]
expr_quintile4 <-comb[comb$Gene %in% bins$Gene[bins$expr_bin=="Bin7" | bins$expr_bin=="Bin8"],]
expr_quintile5 <-comb[comb$Gene %in% bins$Gene[bins$expr_bin=="Bin9" | bins$expr_bin=="Bin10"],]


make_bed_file <- function(geneList) {
  
  toUse <- geneList
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
  return(bed_file)
}

bed_q1 <- make_bed_file(expr_quintile1)
bed_q2 <- make_bed_file(expr_quintile2)
bed_q3 <- make_bed_file(expr_quintile3)
bed_q4 <- make_bed_file(expr_quintile4)
bed_q5 <- make_bed_file(expr_quintile5)

write.table(bed_q5, 
            file = "/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/quintile5_expr_genes_skin.bed", 
            sep="\t", quote = F, row.names = F, col.names = F)

