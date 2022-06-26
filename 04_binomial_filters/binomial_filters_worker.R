## The purpose of this script is to run binomial filters to remove germline variants from somatic mutation calls.

parameters <- commandArgs(TRUE)

patient <- parameters[1]
sex <- parameters[2]
script_dir <- parameters[3]
output_dir <- parameters[4]
pileup_dir <- parameters[5]

#script_dir="/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/04_binomial_filters/"
#script_dir="/Users/so11/phd/psoriasis/scratch/"
#output_dir=paste("/lustre/scratch119/humgen/projects/psoriasis/binomial_filters/", patient, "/", sep="") 
#output_dir="/Users/so11/phd/psoriasis/scratch/"
#pileup_dir <- "/lustre/scratch119/humgen/projects/psoriasis/pileups/"
#pileup_dir <- "/Users/so11/phd/psoriasis/scratch/"
#patient <- "patient05"
#sex <- "male"

source(paste(script_dir, "binomial_filters_function_archive.R", sep=""))

## Hard-coded parameters
gt_matrix_lower <- 0.01     # Everything below this discarded
gt_matrix_upper <- 0.05     # Everything below this set to 'maybe' for tree building
cut_off <- 1                # Sets significant threshold for exact binomial test as q<0.001 (10^-cut_off). Set as 1 for WGS
minimum_mutant_reads <- 0   # Minimum number of supporting reads for calling a mutation. Set as 0 for WGS, 3 for WES
minimum_coverage_of_site <- 0  # Minimum coverage of site to call a mutation at the site. Set as 0 for WGS, 4 for WES. 

options(stringsAsFactors = F)


NR <- read.table(paste(pileup_dir, patient,"/",patient, "_covs.tsv", sep=""), h=T)
NV <- read.table(paste(pileup_dir, patient,"/",patient, "_alts.tsv", sep=""), h=T)

WTR <- NR-NV
VAFs <- NV/NR

## If the coverage at a particular site is 0 then we can get NA in the VAFs matrix
## Set that mutation as 'not present'
VAFs[is.na(VAFs)] <- 0
  
  germline_exact=exact.binomial(NV=NV,NR=NR,gender=sex,cutoff=-cut_off)
  germline_exact[is.na(germline_exact)] <- FALSE
  germline_binom <- beta.binom.filter(NR, NV)
  germline_binom[is.na(germline_binom)] <- FALSE
  
  germline <- !germline_exact & germline_binom
 
  passing_vafs_binom <- VAFs[germline,]
  failing_vafs_binom <- VAFs[!germline,]
  
  NR_pass_binom <- NR[germline,]
  NV_pass_binom <- NV[germline,]
  VAFs_pass_binom <- VAFs[germline,]

  ## I want to next filter based on minimum coverage and the minimum number
  ## of mutant reads. I only want to filter mutations if the minimums are not reached 
  ## in ANY sample. I.e keep low coverage mutations in sample A if I believe they are 
  ## genuine calls based on sample B, which has more data. 

  max_mutant_reads <- apply(NV_pass_binom, 1, max) >= minimum_mutant_reads
  max_coverage <- apply(NR_pass_binom, 1, max) >= minimum_coverage_of_site
  NR_pass <- NR_pass_binom[max_mutant_reads & max_coverage, ]
  NV_pass <- NV_pass_binom[max_mutant_reads & max_coverage, ]
  VAFs_pass <- VAFs_pass_binom[max_mutant_reads & max_coverage, ]
  passing_vafs <- passing_vafs_binom[max_mutant_reads & max_coverage,]
  failing_vafs <- failing_vafs_binom[max_mutant_reads & max_coverage,]


  mutations <- data.frame(chr=unlist(strsplit(rownames(NR_pass), split=":"))[c(T,F,F,F)],
                         pos=as.numeric(unlist(strsplit(rownames(NR_pass), split=":"))[c(F,T,F,F)]),
                         ref=unlist(strsplit(rownames(NR_pass), split=":"))[c(F,F,T,F)],
                         alt=unlist(strsplit(rownames(NR_pass), split=":"))[c(F,F,F,T)])
  
mutations$type <- "indel"
mutations$type[nchar(mutations$ref)==2 & nchar(mutations$alt)==2] <- "dbs"
mutations$type[nchar(mutations$ref)==1 & nchar(mutations$alt)==1] <- "sbs"

  
  write.table(NR_pass, paste(output_dir, patient, "_NR_pass_allMuts.txt", sep=""), quote=F, sep="\t")
  write.table(NV_pass, paste(output_dir, patient, "_NV_pass_allMuts.txt", sep=""), quote=F, sep="\t")
  write.table(VAFs_pass, paste(output_dir, patient, "_VAFs_pass_allMuts.txt", sep=""), quote=F, sep="\t")

  write.table(NR_pass[mutations$type!="indel",], paste(output_dir, patient, "_NR_pass_allSubs.txt", sep=""), quote=F, sep="\t")
  write.table(NV_pass[mutations$type!="indel",], paste(output_dir, patient, "_NV_pass_allSubs.txt", sep=""), quote=F, sep="\t")
  write.table(VAFs_pass[mutations$type!="indel",], paste(output_dir, patient, "_VAFs_pass_allSubs.txt", sep=""), quote=F, sep="\t")

  write.table(NR_pass[mutations$type=="sbs",], paste(output_dir, patient, "_NR_pass_sbs.txt", sep=""), quote=F, sep="\t")
  write.table(NV_pass[mutations$type=="sbs",], paste(output_dir, patient, "_NV_pass_sbs.txt", sep=""), quote=F, sep="\t")
  write.table(VAFs_pass[mutations$type=="sbs",], paste(output_dir, patient, "_VAFs_pass_sbs.txt", sep=""), quote=F, sep="\t")

  write.table(NR_pass[mutations$type=="indel",], paste(output_dir, patient, "_NR_pass_indel.txt", sep=""), quote=F, sep="\t")
  write.table(NV_pass[mutations$type=="indel",], paste(output_dir, patient, "_NV_pass_indel.txt", sep=""), quote=F, sep="\t")
  write.table(VAFs_pass[mutations$type=="indel",], paste(output_dir, patient, "_VAFs_pass_indel.txt", sep=""), quote=F, sep="\t")
  
  genotype_bin <- as.matrix(passing_vafs)
  genotype_bin[genotype_bin<gt_matrix_lower]=0
  genotype_bin[genotype_bin>=gt_matrix_upper]=1
  genotype_bin[genotype_bin>0&genotype_bin<1]=0.5
  
write.table(genotype_bin, file=paste(output_dir, patient, "_genotype_allMuts.txt", sep=""), sep = "\t", quote=F)
write.table(genotype_bin[mutations$type!="indel",], file=paste(output_dir, patient, "_genotype_allSubs.txt", sep=""), sep = "\t", quote=F)
write.table(genotype_bin[mutations$type=="sbs",], file=paste(output_dir, patient, "_genotype_sbs.txt", sep=""), sep = "\t", quote=F)
write.table(genotype_bin[mutations$type=="indel",], file=paste(output_dir, patient, "_genotype_indels.txt", sep=""), sep = "\t", quote=F)
  
  datam <- melt(VAFs_pass)
  datam$value[datam$value<0.01] <- NA
  datam$type <- substr(datam[,1], 4,4)
  datam$type[datam$type=="H"]  <- "Healthy"
  datam$type[datam$type=="L"]  <- "Lesional"

  pdf(paste(output_dir, patient, "_heatmap_and_vafs_allMuts.pdf", sep=""), onefile=T)
    plot <- heatmap(genotype_bin,scale='none',
          col=c("aliceblue","lightblue","steelblue"),mar=c(8,8))
    print(plot)

    plot <- ggplot(datam, aes(x=value, fill=type)) + geom_histogram() + facet_wrap(~variable, nrow=3) + 
    labs(x="VAF", y="Number of mutations") + scale_fill_manual(values=c("#FFDBAC","#e18377")) 
    print(plot)
  dev.off()
  

  
