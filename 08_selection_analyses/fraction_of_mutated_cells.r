
## The purpose of this script is to determine the fraction of mutations of mutated
## cells for each individual and overall. 
library(reshape2)
library(ggplot2)

genes2test <- c("NOTCH1", "NOTCH2", "PPM1D","FAT1", "TP53", "ZFP36L2", 
                "MRC1", "GXYLT1", "CHEK2", "EEF1A1")

all_muts <- read.table("/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/annotated_coding_mutations.txt", h=T, stringsAsFactors = F)

patients <- substr(all_muts$sampleID, 2,5)
patients <- gsub("H", "", patients)
patients <- gsub("L", "", patients)
patients <- gsub("_", "", patients)
all_muts$Patient <- paste("patient", patients, sep="")

sample_meta <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/sample_info/sample_meta.txt", h=T, stringsAsFactors = F)
patient_meta <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/sample_info/patient_meta.txt", h=T, stringsAsFactors = F)
patients <- unique(sample_meta$patient_ID[sample_meta$Seq_status=="WES" & is.na(sample_meta$excl_criteria)])


muts_of_interest <- all_muts[all_muts$gene %in% genes2test & all_muts$impact!="Synonymous",]
muts_of_interest$mutID <- paste(muts_of_interest$chr, muts_of_interest$pos, muts_of_interest$ref, muts_of_interest$mut, sep=":")


fraction_mutated <- matrix(data=NA, nrow=length(patients), ncol=length(genes2test))
i <- 1
for(patient in patients) {
  VAF_matrix <- read.table(paste("/lustre/scratch119/humgen/projects/psoriasis/binomial_filters/", patient, "/",patient,"_VAFs_pass_allMuts.txt", sep=""), h=T)
  rownames(VAF_matrix) <- gsub("chr", "", rownames(VAF_matrix))
  j <- 1
  for(gene in genes2test) {
    muts_in_gene <- muts_of_interest$mutID[muts_of_interest$Patient==patient & muts_of_interest$gene==gene]
    if(length(muts_in_gene)==0) {
      # Patient has no mutations in the gene. Fraction of mutated cells is zero.
      fraction_mutated[i,j] <- 0
    } else {
      # For each mutation, the fraction of cells carrying the mutation is the sum of two times the VAF of the mutation in each sample.
      # For now I'm going to assume that each clone only carries a single mutation in each gene so the total number of mutated cells
      # is the sum of 2*VAFs*volume of each sample divided by the total sample volume
      # This assumption actually doesn't hold true. Some samples have multiple mutations in the same genes. See for example patient17 and FAT1.
      
      ## Note: If any genes on the sex-chromosomes are under positive selection, that will have to be accounted for separately. 
      ## Note: In time, only include samples where the coverage of the gene in question is such that the coverage*mVAF > 4. Those are the samples for which 
      ## we would be able to call mutation in that gene. 
      x <- VAF_matrix*sample_meta$Volume[sample_meta$sampleID %in% colnames(VAF_matrix)]
      fraction_mutated[i,j] <- sum(rowSums(VAF_matrix[muts_in_gene,])*2)/ncol(VAF_matrix)  ## without weighing by sample volume
      #fraction_mutated[i,j] <- sum(rowSums(x[muts_in_gene,])*2)/sum(sample_meta$Volume[sample_meta$sampleID %in% colnames(VAF_matrix)])
    }
    j <- j+1
  }
  i <- i+1
}

rownames(fraction_mutated) <- patients
colnames(fraction_mutated) <- genes2test
fraction_mutated_df <- data.frame(fraction_mutated)
fraction_mutated_df$Patient <- rownames(fraction_mutated_df)
patient_meta <- merge(fraction_mutated_df, patient_meta, by.x="Patient", by.y="Patient.ID")

patient_meta$Age_at_sampling <- patient_meta$Year_of_sampling-patient_meta$Year_of_birth
patient_meta$Total_driver_prevalence <- patient_meta$NOTCH1+patient_meta$NOTCH2+patient_meta$PPM1D+patient_meta$TP53+patient_meta$FAT1 + patient_meta$TP53
patient_meta$PUVA_seen <- "No"
patient_meta$PUVA_seen[patient_meta$Patient %in% c("patient14","patient18","patient21","patient24","patient27","patient06","patient34")] <- "Yes"
patient_meta$Disease_duration <- patient_meta$Year_of_sampling - patient_meta$Year_of_onset

x <- melt(fraction_mutated_df)
x <- merge(x, patient_meta[, c("Patient", "Age_at_sampling")])

p1 <- ggplot(x, aes(x=reorder(Patient, Age_at_sampling), y=variable)) + theme_classic(base_size = 15) +
  geom_tile(aes(fill=value)) + theme(axis.text.x = element_text(angle=90)) + labs(x="", y="") + 
  #geom_text(aes(label = round(value,2))) + 
  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "none")


p2 <- ggplot(patient_meta, aes(x=reorder(Patient, Age_at_sampling), y="")) + theme_void() +
  geom_tile(aes(fill=Sex)) + labs(x="", y="") + 
  scale_fill_manual(values=c("pink","#63D2FF")) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
  theme(legend.position = "none")

p3 <- ggplot(patient_meta, aes(x=reorder(Patient, Age_at_sampling), y="")) + theme_void() +
  geom_tile(aes(fill=as.factor(Current_therapy._topical_steroids))) + theme(axis.text.x = element_text(angle=90)) + labs(x="", y="") + 
  scale_fill_manual(values=c("#27187E","#F79824")) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none")

p4 <- ggplot(patient_meta, aes(x=reorder(Patient, Age_at_sampling), y="")) + theme_void() +
  geom_tile(aes(fill=PUVA_seen)) + theme(axis.text.x = element_text(angle=90)) + labs(x="", y="") + 
  scale_fill_manual(values=c("#AFBFC0","#BA5624")) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none")

p5 <- ggplot(patient_meta, aes(x=reorder(Patient, Age_at_sampling), y="")) + theme_void() +
  geom_tile(aes(fill=Disease_duration)) + theme(axis.text.x = element_text(angle=90)) + labs(x="", y="") + 
  scale_fill_gradient(low = "white", high = "#FDCA40") +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none")

p6 <- ggplot(patient_meta, aes(x=reorder(Patient, Age_at_sampling), y="")) + theme_void() +
  geom_tile(aes(fill=Age_at_sampling)) + theme(axis.text.x = element_text(angle=90)) + labs(x="", y="") + 
  scale_fill_gradient(low = "white", high = "#ED254E") +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none")

plot_grid( p4, p3,p2,p5,p6,p1, rel_heights = c(1,1,1,1,1,10), nrow=6, align="v")


df <- melt(fraction_mutated)

f1 <- ggplot(df[df$value<1,], aes(x=value*100, fill=Var2)) + geom_histogram() + 
  facet_wrap(.~Var2, nrow=2) + theme_classic(base_size = 20) +
  labs(x="Percentage of mutated cells", y="Number of patients") +
  scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")

f2 <- ggplot(patient_meta[patient_meta$Patient!="patient17",], aes(x=Age_at_sampling, y=Total_driver_prevalence)) +
  geom_point(size=3) + theme_classic(base_size = 20) + labs(x="Age", y="Sum of driver fractions \n accross genes") +
  geom_smooth(method="lm") + ylim(c(0,1))


###
##
##  COMPARE THE FRACTION OF MUTATED CELLS BETWEEN LESIONAL AND NON-LESIONAL SKIN
##
################################################################################


fraction_mutated_L <- matrix(data=NA, nrow=length(patients), ncol=length(genes2test))
fraction_mutated_NL <- matrix(data=NA, nrow=length(patients), ncol=length(genes2test))
nr_lesional_samples <- as.numeric()
nr_healthy_samples <- as.numeric()
i <- 1
for(patient in patients) {
  VAF_matrix <- read.table(paste("/lustre/scratch119/humgen/projects/psoriasis/binomial_filters/", patient, "/",patient,"_VAFs_pass_allMuts.txt", sep=""), h=T)
  
  j <- 1
  for(gene in genes2test) {
    muts_in_gene <- muts_of_interest$mutID[muts_of_interest$Patient==patient & muts_of_interest$gene==gene]
    if(length(muts_in_gene)==0) {
      # Patient has no mutations in the gene. Fraction of mutated cells is zero.
      fraction_mutated_L[i,j] <- 0
      fraction_mutated_NL[i,j] <- 0
    } else {
      # For each mutation, the fraction of cells carrying the mutation is the sum of two times the VAF of the mutation in each sample.
      # For now I'm going to assume that each clone only carries a single mutation in each gene so the total number of mutated cells
      # is the sum of 2*VAFs*volume of each sample divided by the total sample volume
      # This assumption actually doesn't hold true. Some samples have multiple mutations in the same genes. See for example patient17 and FAT1.
      
      ## Note: If any genes on the sex-chromosomes are under positive selection, that will have to be accounted for separately. 
      ## Note: In time, only include samples where the coverage of the gene in question is such that the coverage*mVAF > 4. Those are the samples for which 
      ## we would be able to call mutation in that gene. 
      lesional_samples <- colnames(VAF_matrix) %in% sample_meta$sampleID[sample_meta$type=="Lesional"]
      healthy_samples <- colnames(VAF_matrix) %in% sample_meta$sampleID[sample_meta$type %in% c("Non-lesional", "Non_lesional", "Healthy")]
      
      if(sum(lesional_samples)>1) {
        fraction_mutated_L[i,j] <- sum(rowSums(VAF_matrix[muts_in_gene,lesional_samples])*2)/sum(colnames(VAF_matrix) %in% sample_meta$sampleID[sample_meta$type=="Lesional"])
      }
      if(sum(healthy_samples)>1) {
        fraction_mutated_NL[i,j] <- sum(rowSums(VAF_matrix[muts_in_gene,healthy_samples])*2)/sum(colnames(VAF_matrix) %in% sample_meta$sampleID[sample_meta$type %in% c("Non-lesional", "Non_lesional", "Healthy")])
      } else if(sum(healthy_samples)==1) {
        fraction_mutated_NL[i,j] <- sum(sum(VAF_matrix[muts_in_gene,healthy_samples])*2)/sum(colnames(VAF_matrix) %in% sample_meta$sampleID[sample_meta$type %in% c("Non-lesional", "Non_lesional", "Healthy")])
      } else {
        fraction_mutated_NL[i,j] <- NA
      }
    }
    j <- j+1
  }
  nr_healthy_samples[i] <-sum(healthy_samples)
  nr_lesional_samples[i] <- sum(lesional_samples)
  i <- i+1
}

rownames(fraction_mutated_L) = rownames(fraction_mutated_NL) <- patients
colnames(fraction_mutated_L) = colnames(fraction_mutated_NL) <- genes2test

fraction_mutated_L <- data.frame(fraction_mutated_L)
fraction_mutated_NL <- data.frame(fraction_mutated_NL)
fraction_mutated_L$type <- "Lesional"
fraction_mutated_NL$type <- "Non-Lesional"

test <- melt(rbind(fraction_mutated_L, fraction_mutated_NL))
f3 <- ggplot(test, aes(x=variable, y=value, fill=type)) + geom_boxplot() + ylim(c(0,1)) + 
  theme_classic(base_size = 20) + scale_fill_manual(values=c(col_lesion,col_nonLes)) +
  theme(axis.text.x = element_text(angle=90)) + labs(y="Fraction of mutated cells", x="") + 
  theme(legend.title = element_blank(), legend.position = "top")

plot_grid(f1,plot_grid(f3,f2, nrow=2), ncol=2)


wilcox.test(fraction_mutated_L$NOTCH1, fraction_mutated_NL$NOTCH1, paired = TRUE)
wilcox.test(fraction_mutated_L$NOTCH2, fraction_mutated_NL$NOTCH2, paired = TRUE)
wilcox.test(fraction_mutated_L$PPM1D, fraction_mutated_NL$PPM1D, paired = TRUE)
wilcox.test(fraction_mutated_L$FAT1, fraction_mutated_NL$FAT1, paired = TRUE)
wilcox.test(fraction_mutated_L$TP53, fraction_mutated_NL$TP53, paired = TRUE)
wilcox.test(fraction_mutated_L$ZFP36L2, fraction_mutated_NL$ZFP36L2, paired = TRUE)


wilcox.test(fraction_mutated_L$NOTCH1[!(fraction_mutated_L$NOTCH1==0 & fraction_mutated_NL$NOTCH1==0)], 
            fraction_mutated_NL$NOTCH1[!(fraction_mutated_L$NOTCH1==0 & fraction_mutated_NL$NOTCH1==0)], paired = TRUE)
wilcox.test(fraction_mutated_L$NOTCH2[!(fraction_mutated_L$NOTCH2==0 & fraction_mutated_NL$NOTCH2==0)], 
            fraction_mutated_NL$NOTCH2[!(fraction_mutated_L$NOTCH2==0 & fraction_mutated_NL$NOTCH2==0)], paired = TRUE)
wilcox.test(fraction_mutated_L$PPM1D[!(fraction_mutated_L$PPM1D==0 & fraction_mutated_NL$PPM1D==0)], 
            fraction_mutated_NL$PPM1D[!(fraction_mutated_L$PPM1D==0 & fraction_mutated_NL$PPM1D==0)], paired = TRUE)
wilcox.test(fraction_mutated_L$FAT1[!(fraction_mutated_L$FAT1==0 & fraction_mutated_NL$FAT1==0)], 
            fraction_mutated_NL$FAT1[!(fraction_mutated_L$FAT1==0 & fraction_mutated_NL$FAT1==0)], paired = TRUE)
wilcox.test(fraction_mutated_L$TP53[!(fraction_mutated_L$TP53==0 & fraction_mutated_NL$TP53==0)], 
            fraction_mutated_NL$TP53[!(fraction_mutated_L$TP53==0 & fraction_mutated_NL$TP53==0)], paired = TRUE)
wilcox.test(fraction_mutated_L$ZFP36L2[!(fraction_mutated_L$ZFP36L2==0 & fraction_mutated_NL$ZFP36L2==0)], 
            fraction_mutated_NL$ZFP36L2[!(fraction_mutated_L$ZFP36L2==0 & fraction_mutated_NL$ZFP36L2==0)], paired = TRUE)


wilcox.test(rowSums(fraction_mutated_L[,c(1:6)]), rowSums(fraction_mutated_NL[,c(1:6)]), paired = TRUE)

