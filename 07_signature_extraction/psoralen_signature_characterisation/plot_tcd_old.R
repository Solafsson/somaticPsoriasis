
## The purpose of this script is to read in the output from 
## puva_tcd.sh and plot the results. 


Tmut_transcribed <- read.table("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/combined_Tmut_transcribed.bed")
Amut_transcribed <- read.table("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/combined_Amut_transcribed.bed")
#Cuv_transcribed <- read.table("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/combined_Cuv_transcribed_sbs7b.bed")
#Guv_transcribed <- read.table("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/combined_Guv_transcribed_sbs7b.bed")
#Cuv_transcribed <- read.table("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/wes_Cuv_transcribed.bed")
#Guv_transcribed <- read.table("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/wes_Guv_transcribed.bed")
#Cother_transcribed <- read.table("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/combined_Cother_transcribed.bed")
#Gother_transcribed <- read.table("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/combined_Gother_transcribed.bed")



Tmut_transcribed$V6 <- factor(Tmut_transcribed$V6, levels=c("-10kb", "-9kb", "-8kb", "-7kb", "-6kb", "-5kb", "-4kb", "-3kb", "-2kb", "-1kb",
                                          "+1kb", "+2kb", "+3kb", "+4kb", "+5kb", "+6kb", "+7kb", "+8kb", "+9kb", "+10kb"))
Amut_transcribed$V6 <- factor(Amut_transcribed$V6, levels=c("-10kb", "-9kb", "-8kb", "-7kb", "-6kb", "-5kb", "-4kb", "-3kb", "-2kb", "-1kb",
                                            "+1kb", "+2kb", "+3kb", "+4kb", "+5kb", "+6kb", "+7kb", "+8kb", "+9kb", "+10kb"))


Tmut_transcribed$type <- "T>[ACG] transcribed, A>[CGT] non-transcribed"
Amut_transcribed$type <- "A>[CGT] transcribed, T>[ACG] non-transcribed"


comb <- rbind(Tmut_transcribed, Amut_transcribed,Cuv_transcribed,Guv_transcribed)

library(ggplot2)
library(reshape2)

m <- melt(data.frame(table(comb$V6, comb$type)))
m$variable <- NULL
colnames(m) <- c("DistanceFromTSS", "MutationType", "MutationCount")

standardized <- as.numeric()
for(i in 1:nrow(m)) {
  denom <- m$MutationCount[m$DistanceFromTSS=="-10kb" & m$MutationType==m$MutationType[i]]
  standardized[i] <- m$MutationCount[i]/denom
}
m$standardized <- standardized

write.table(m, "/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/puva_tcd.txt", row.names = F, quote=F, sep="\t")

ggplot(m, aes(x=DistanceFromTSS, y=standardized, colour=MutationType)) + geom_point() + theme_bw(base_size = 14) +
  geom_line(aes(x=DistanceFromTSS,y=standardized, color = MutationType, group = MutationType)) + 
  labs(x="Distance from the TSS", y="Mutation rate relative to intergenic regions", colour="") + 
  theme(legend.position = "top")


## Test a step function. Want to see if the upwards trend of the A>[CGT] Transcribed
## T>[ACG] non-transcribed is significantly different from a straight line. 
counts <- data.frame(table(Amut_transcribed$V6))
counts$upstream_TSS <- 1
counts$upstream_TSS[grep("-", as.character(counts$Var1))] <- 0

null_model <- lm(Freq ~ 1, data=counts)
alt_model <- lm(Freq ~ I(counts$upstream_TSS), data=counts)

anova(null_model, alt_model)$"Pr(>F)"[2]







