
## The purpose of this script is to read in the output from 
## puva_tcd.sh and plot the results. 

for(q in c(1:5)) {
  assign(paste("quintile", q, "_Tmut_transcribed", sep=""), 
         read.table(paste("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/quintile", q, "_Tmut_transc.bed", sep="")))

  assign(paste("quintile", q, "_Amut_transcribed", sep=""), 
         read.table(paste("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/quintile", q, "_Amut_transc.bed", sep="")))

}

quintile1_Amut_transcribed$quintile=quintile1_Tmut_transcribed$quintile <- "Expression quintile1"
quintile2_Amut_transcribed$quintile=quintile2_Tmut_transcribed$quintile <- "Expression quintile2"
quintile3_Amut_transcribed$quintile=quintile3_Tmut_transcribed$quintile <- "Expression quintile3"
quintile4_Amut_transcribed$quintile=quintile4_Tmut_transcribed$quintile <- "Expression quintile4"
quintile5_Amut_transcribed$quintile=quintile5_Tmut_transcribed$quintile <- "Expression quintile5"

quintile1_Amut_transcribed$type=quintile2_Amut_transcribed$type=quintile3_Amut_transcribed$type=quintile4_Amut_transcribed$type=quintile5_Amut_transcribed$type <- "A>[CGT] transcribed,\nT>[ACG] non-transcribed" 
quintile1_Tmut_transcribed$type=quintile2_Tmut_transcribed$type=quintile3_Tmut_transcribed$type=quintile4_Tmut_transcribed$type=quintile5_Tmut_transcribed$type <- "T>[ACG] transcribed,\nA>[CGT] non-transcribed" 

comb <- rbind(quintile1_Amut_transcribed,quintile1_Tmut_transcribed,quintile2_Amut_transcribed,quintile2_Tmut_transcribed,
              quintile3_Amut_transcribed,quintile3_Tmut_transcribed,quintile4_Amut_transcribed,
              quintile4_Tmut_transcribed,quintile5_Amut_transcribed,quintile5_Tmut_transcribed)

comb$V6 <- factor(comb$V6, levels=c("-10kb", "-9kb", "-8kb", "-7kb", "-6kb", "-5kb", "-4kb", "-3kb", "-2kb", "-1kb",
                                    "+1kb", "+2kb", "+3kb", "+4kb", "+5kb", "+6kb", "+7kb", "+8kb", "+9kb", "+10kb"))

m <- melt(data.frame(table(comb$V6, comb$quintile, comb$type)))
m$variable <- NULL
colnames(m) <- c("DistanceFromTSS","expression", "MutationType", "MutationCount")

standardized <- as.numeric()
for(i in 1:nrow(m)) {
  denom <- m$MutationCount[m$DistanceFromTSS=="-10kb" & m$MutationType==m$MutationType[i] & m$expression==m$expression[i]]
  standardized[i] <- m$MutationCount[i]/denom
}
m$standardized <- standardized

m$class <- paste(m$expression, m$MutationType)

write.table(m, "/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/puva_tcd.txt", row.names = F, quote=F, sep="\t")

ggplot(m, aes(x=DistanceFromTSS, y=standardized, colour=expression)) + geom_point() + theme_bw(base_size = 14) +
  geom_line(aes(x=DistanceFromTSS,y=standardized, color = expression, group = class, linetype=MutationType)) + 
  labs(x="Distance from the TSS", y="Mutation rate relative to intergenic regions", colour="", linetype="") + 
  scale_y_continuous(trans="log2") + scale_colour_manual(values=c("#F79256","#FBD1A2", "#7DCFB6","#00B2CA","#1D4E89"))
  




