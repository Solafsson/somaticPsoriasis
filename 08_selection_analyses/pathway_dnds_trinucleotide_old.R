## Read in a very slightly modified version of the dndscv() function which lumps
## Nonsense and splice mutations together under "truncating" mutations. 
source("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/08_selection_analyses/dndscv_pvalpois.r")

## These gene lists were defined a priory, without us having looked at the data first. 
## See the methods section of the manuscript for details of how pathway membership was
## determined. 
geneLists <- c("Normal_skin_pos", "BCC","GWAS_psoriasis","IL17","IL12_23","TNF","IL36_MyD88","IFNg","MHC_classI","TLR", "IBD_mucosa")
for(geneL in geneLists) {
  assign(paste(geneL, "_genes", sep=""), read.table(paste("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/08_selection_analyses/pathway_dnds_gene_lists/", geneL, ".txt", sep="")))
}

## This analysis doesn't include indels/DBS mutations so there is no need to split the
## mutation list like I did above. It doesn't matter which df is used. 
for(geneL in geneLists) {
  assign(paste("dnds_",geneL,sep=""), dndscv_pvalspois(dbs_only, get(paste(geneL, "_genes", sep=""))$V1,max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf,refdb=refcds_38, cv=scores))
}

df <- data.frame()
for(geneL in geneLists) {
  assign(paste("sel_cv_", geneL, sep=""), get(paste("dnds_",geneL,sep=""))$sel_cv)
  
  df_tmp <- data.frame(get(paste("dnds_",geneL,sep=""))$globaldnds)
  df_tmp$geneList <- geneL
  df <- rbind(df, df_tmp)
}
mis_and_trunc <- df[df$name %in% c("wmis", "wtru"),]

P_mis <- as.numeric()
P_trunc <- as.numeric()

for(geneL in geneLists) {
  P_mis <- c(P_mis, summary(get(paste("dnds_", geneL, sep=""))$poissmodel)$coefficients[c("wmis"),4] )
  P_trunc <- c(P_trunc, summary(get(paste("dnds_", geneL, sep=""))$poissmodel)$coefficients[c("wtru"),4] )
}
mis_and_trunc$p_pois <- NA
mis_and_trunc$p_pois[seq(from=1, to=nrow(mis_and_trunc), by=2)] <- P_mis
mis_and_trunc$p_pois[seq(from=2, to=nrow(mis_and_trunc), by=2)] <- P_trunc
mis_and_trunc$q_BHcorr <- p.adjust(mis_and_trunc$p_pois, method="BH")

mis_and_trunc
mis_and_trunc[order(mis_and_trunc$q_BHcorr),]
#write.table(mis_and_trunc, "/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/pathway_dnds_mis_and_trunc_w_covs.txt", row.names=F, quote=F, sep="\t")

