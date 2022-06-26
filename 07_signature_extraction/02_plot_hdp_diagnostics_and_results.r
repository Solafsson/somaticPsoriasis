

hdp_posterior_dir="/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/hdp/snv/"
output_dir <- "/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/hdp/snv/"
# read in PCAWG sigs
all_sigs <- read.csv("/lustre/scratch119/humgen/projects/psoriasis/resources/sigProfiler_SBS_signatures_summing_to_one.csv", header = T, row.names = 1, stringsAsFactors = F)

#.libPaths("/lustre/scratch114/projects/crohns/somatic_ibd_p1/R_packages_farm5_R3.6_install/")
.libPaths("/lustre/scratch119/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")
library(hdp)

gdsigs <- c("SBS1", "SBS5","SBS2", "SBS13","SBS7a","SBS7b","SBS7c","SBS7d", "SBS17a", "SBS17b","SBS18","SBS38")

load(paste(output_dir, "branch_mutational_matrix.txt", sep=""))
tnt <- mut_mat_t
newColNames <- paste(substr(colnames(mut_mat_t), 3,3), ".", substr(colnames(mut_mat_t), 5,5), ".in.", substr(colnames(mut_mat_t), 1,1), substr(colnames(mut_mat_t), 3,3), substr(colnames(mut_mat_t), 7,7), sep="")
colnames(tnt) <- newColNames
patients <- sapply(strsplit(rownames(tnt), "_"), "[[", 1)



worked <- c(1:20)

for(i in worked) {
  load(paste(hdp_posterior_dir, "hdpout_", i, ".RData", sep=""))
}

chlist <- vector("list", length(worked))
for(i in 1:length(worked)) {
  chlist[[i]] <- get(paste("hdp_", worked[i], sep=""))
}


luad_multi <- hdp_multi_chain(chlist)

# plot diagnostics
pdf(paste(output_dir,"hdp_diagnostics.pdf", sep=""), onefile=T)
    par(mfrow=c(4,5), mar=c(5,4,1,1))
    lapply(chains(luad_multi), plot_lik, bty='L', start=1000)
    lapply(chains(luad_multi), plot_lik, bty='L')
    lapply(chains(luad_multi), plot_numcluster, bty='L')
    lapply(chains(luad_multi), plot_data_assigned, bty='L')
dev.off()




## extract consensus components / signatures
luad_multi <- hdp_extract_components(luad_multi, cos.merge = 0.9)
luad_multi
numcomp(luad_multi)
prop.ex(luad_multi)
save(luad_multi, file="/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/wes_hdp_components.RData")

par(mfrow=c(1,1), mar=c(4, 4, 2, 2))

# plot number of mutations per signature (one dot per posterior sample)
pdf(paste(output_dir,"hdp_nr_muts_per_signature.pdf", sep=""), onefile=T)
    plot_comp_size(luad_multi, bty="L")
dev.off()

## Get the counts and the fraction on a numeric format:
sums <- t(sapply(comp_categ_counts(luad_multi), rowSums))
totalMuts <- sum(apply(sums,1,mean))
apply(sums,1,mean)/totalMuts * 100


# plot components / signatures
# pick your colours
#mut_colours <- c(RColorBrewer::brewer.pal(10, 'Paired')[seq(1,10,2)], 'grey70')
## Use same colours as COSMIC
mut_colours <- c("#16bdebff", "#000000ff", "#e22926ff", "#999999ff", "#9fce67ff", "#ecc6c5ff")
# labels along bototm x-axis
trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
# group labels along the top (and controls colour grouping)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))

par(mfrow=c(5,2), mar=c(4, 4, 2, 2))
## 


pdf(paste(output_dir,"hdp_components.pdf", sep=""), onefile=T)
plot_comp_distn(luad_multi, cat_names=trinuc_context,
                grouping=group_factor, col=mut_colours,
                col_nonsig="grey80", show_group_labels=TRUE,incl_comp0=F,
                cred_int=F,plot_title=NULL)
dev.off()

## Once you've figured out which signature is which then you can
## make the plot with signature names, but this will have to be done
## manually as the number and order of signatures will change.

pdf(paste(output_dir,"hdp_components_signatures.pdf", sep=""), onefile=T)
  plot_comp_distn(luad_multi, cat_names=trinuc_context,
                  grouping=group_factor, col=mut_colours,
                   show_group_labels=TRUE,
                  plot_title = c("Unassigned","SBS7b - UV","PUVA treatment", "SBS1/5", "UV-component N1",
                                 "Unknown-N2","SBS13 - APOBEC", "Unknown-N3","SBS7d - UV","SBS2 - APOBEC"
                                 ), cred_int=F)

dev.off()


png(paste("/lustre/scratch119/humgen/projects/psoriasis/scratch/","hdp_components_signatures.png", sep=""), width = 1300, height=1100)
par(mfrow=c(5,2), mar=c(4, 4, 2, 2))
plot_comp_distn(luad_multi, cat_names=trinuc_context,
                grouping=group_factor, col=mut_colours,
                show_group_labels=TRUE,
                plot_title = c("Unassigned","SBS7b - UV","PUVA treatment", "SBS1/5", "UV-component N1",
                               "UV-component N2","SBS2 - APOBEC","SBS7d - UV", "New component N3",
                               "SBS13 - APOBEC"), cred_int=F)

dev.off()


par(mfrow=c(1,1), mar=c(4, 4, 2, 2))
## The object contains 1 root node, one node that is the parent of all data, 1 pseuodo node for each prior
## and one node for each patient. 
## 1+1+length(gdsigs) + length(unique(patients))
## The first node containing sample information then is the next node after that. 
## This number will change as more samples and a different number of priors are used so beware. 
first_sample_node <- 1+1+length(gdsigs) + length(unique(patients)) + 1
last_sample_node <- length(comp_dp_distn(luad_multi)$cred.int)

## Need to come up with a better solution for the color palette. Need a palette that can handle more colours. 
pdf(paste(output_dir,"hdp_signature_exposure.pdf", sep=""), onefile=T)
plot_dp_comp_exposure(luad_multi, main_text="Signatures - All patients",
                      dpindices=first_sample_node:last_sample_node,
                      col=c(RColorBrewer::brewer.pal(12, "Set3"), "navy", "forestgreen", "red", "yellow", "pink", "brown"),
                      incl_nonsig=FALSE,
                      ylab_numdata = 'SNV count', ylab_exp = 'Signature exposure',
                      leg.title = 'Signature',incl_numdata_plot = T, dpnames = rownames(tnt),
                      las=2,mar=c(3,4,2,0.5))
dev.off()


# To extract the probability distribution for a given signature
comp_distn <- comp_categ_distn(luad_multi)
sig <- comp_distn$mean["N1", ]
names(sig) <- trinuc_context



# to extract relative contributions
toutmeans <- (comp_dp_distn(luad_multi)$mean)
toutmeans <- toutmeans[first_sample_node:(nrow(toutmeans)),]
rownames(toutmeans) <- rownames(tnt)
colnames(toutmeans) <- c("Unassigned", "SBS7b", "PUVA", "SBS1/5","Unknown-N1", "Unknown-N2","SBS2", "Unknown-N3","SBS7c","SBS13")

library(corrplot)
m <- cor(toutmeans)
corrplot(m, method="circle", insig="blank")

write.table(data.frame(toutmeans),paste(output_dir, "branch_exposures_w_prior.txt", sep=""), row.names = T, quote=F)

