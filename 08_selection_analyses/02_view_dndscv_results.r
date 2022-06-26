
mutation_list <- read.table("/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/annotated_coding_mutations.txt", h=T)
sel_cv <- read.table("/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/sel_cv_indel_test.txt", h=T)

signif_genes = sel_cv[sel_cv$qglobal_cv<0.1 | sel_cv$qallsubs_cv<0.05, c(1:6, 17,19)]
mutations_in_signif_genes <- mutation_list[mutation_list$gene %in% signif_genes$gene_name,]

sel_cv_noCov <- read.table("/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/sel_cv_nocov.txt", h=T)
sel_cv_noHype <- read.table("/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/sel_cv_noHype.txt", h=T)
sel_cv_lesional <- read.table("/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/sel_cv_lesional.txt", h=T)

for(gene in unique(mutations_in_signif_genes$gene)) {
  thisGene <- mutations_in_signif_genes[mutations_in_signif_genes$gene==gene, ]
  
  thisGene$Cancer_Type <- "Skin"
  thisGene$End_position <- thisGene$pos
  for(i in 1:nrow(thisGene)) {
    longest <- max(nchar(thisGene$ref[i]), nchar(thisGene$mut[i]))
    thisGene$End_position[i] <- thisGene$pos[i] + longest -1
  }
  
  BioPortal_input <- thisGene[,c("sampleID", "Cancer_Type", "chr", "pos", "End_position", "ref", "mut")]
  colnames(BioPortal_input) <- c("Sample_ID","Cancer_Type","Chromosome","Start_Position","End_Position","Reference_Allele","Variant_Allele")
  
  write.table(BioPortal_input, file=paste("/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/", gene, "_BioPortal_input.txt", sep=""),row.names = F, quote = F, sep="\t")
  
}


mis_and_trunc <- read.table("/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/pathway_dnds_mis_and_trunc_w_covs.txt", h=T)


library(reshape2)
library(ggplot2)
library(scales)
## Get rid of 0 and INF for IL6 (No mutations found)
mis_and_trunc$cilow[mis_and_trunc$cilow<0.1] <- min(mis_and_trunc$cilow[mis_and_trunc$cilow>0.1])
mis_and_trunc$cihigh[mis_and_trunc$cihigh>1000] <- max(mis_and_trunc$cihigh[mis_and_trunc$cihigh<1000])
mis_and_trunc$mle[mis_and_trunc$mle<0.1] <- 1

mis_and_trunc$name <- factor(mis_and_trunc$name, levels=c("wmis", "wtru"), labels = c("Missense", "Truncating"))
mis_and_trunc$geneList <- factor(mis_and_trunc$geneList, levels=c("Normal_skin_pos", "BCC", "GWAS_psoriasis","IBD_mucosa", "TNF", "IFNg", "IL12_23", "IL36_MyD88",
                                                                  "TLR","IL17", "MHC_classI"),
                                 c("Normal skin/SCCs", "BCCs", "Psoriasis GWAS hits", "IBD mucosa", "TNF", "IFNg", "IL12/23", "IL36/MyD88", "TLR", "IL17", "MHC-class I"))




## Make the plot without missense variants 
trunc <- subset(mis_and_trunc, name!="Missense")

dnds_plot_trunc <- ggplot(trunc, aes(x = geneList, y = mle, ymin = cilow, ymax = cihigh)) + 
  geom_point(size=4, position=position_dodge(width=0.5)) + 
  geom_errorbar(size=1, position = "dodge", width=0.5) +theme_classic(base_size = 14) + 
  geom_hline(yintercept = 1, aes(size=2)) + theme(axis.text.x = element_text(angle=75, hjust=1)) +labs(x="", y="dN/dS") + 
  scale_y_continuous(trans=log2_trans()) 

mis <- subset(mis_and_trunc, name=="Missense")
dnds_plot_mis <- ggplot(mis, aes(x = geneList, y = mle, ymin = cilow, ymax = cihigh)) + 
  geom_point(size=4, position=position_dodge(width=0.5)) + 
  geom_errorbar(size=1, position = "dodge", width=0.5) +theme_classic(base_size = 14) + 
  geom_hline(yintercept = 1, aes(size=2)) + theme(axis.text.x = element_text(angle=75, hjust=1)) +labs(x="", y="dN/dS") + 
  scale_y_continuous(trans=log2_trans()) 

ggplot(mis_and_trunc, aes(x=geneList, y=mle, ymin = cilow, ymax = cihigh, colour=name)) +
  geom_point(size=4, position=position_dodge(width=0.5)) + 
  geom_errorbar(size=1, position = "dodge", width=0.5) +theme_classic(base_size = 14) + 
  geom_hline(yintercept = 1, aes(size=2)) + theme(axis.text.x = element_text(angle=75, hjust=1)) +labs(x="", y="dN/dS") + 
  scale_y_continuous(trans=log2_trans()) + theme(legend.position = "top") + labs(colour="")
  