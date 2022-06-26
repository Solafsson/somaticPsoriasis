
## This script contains a bunch of code which was originally in Run_Dirichlet_clustering_posthoc.R
## It depends on the output of that script being available. 

#############################################
## -- Receive command line args --
args = commandArgs(trailingOnly=TRUE)
args

repo_location = args[1]
output_dir = args[2]
patient = args[3]

#############################################
## -- Library loading and sourcing --

.libPaths("/lustre/scratch119/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")

library(RColorBrewer)
library(label.switching)
library(philentropy)
library(ggplot2)
library(GGally)
library(dplyr)
library(data.table)
library(here)
library(farver)

# source(file.path(repo_location, 'lib', 'Mutation_Cluster_Dirichlet_Gibbs_sampler_functions.R'))
# source(file.path(repo_location, 'lib', 'Mutation_Cluster_Dirichlet_posthoc_merge_fx.R'))
source(file.path(repo_location,'Mutation_Cluster_Dirichlet_Gibbs_sampler_functions.R'))
source(file.path(repo_location,'Mutation_Cluster_Dirichlet_posthoc_merge_fx.R'))



### Read in the output from Run_Dirichlet_clustering_posthoc.R
load(paste(output_dir,"ndp_",patient, "_2022_02_15/Rsession_ls.dat", sep=""))

par_min.threshold = 15

#############################################

# Merge clusters which differ only at level of sequencing errors

lymph.merge.ls <- merge.clusters.differ.by.seq.error(lymph.gs, lymph.ls, min.true.VAF = 0.025, overlap.frac=0.025, min.num.muts = par_min.threshold, ls.type="ECR", number.burn.in, version=2)
ls_tbl = tibble(clust_assign = as.numeric(lymph.merge.ls$final.cluster)) %>%
  mutate(mut_id = row_number())
ls_tbl$chrom = lymph.depth$chrom
ls_tbl$pos = lymph.depth$pos
fwrite(ls_tbl, file.path(output_dir, 'clust_assign_posthoc.csv'))
saveRDS(lymph.merge.ls, file.path(output_dir, 'object_lymph.merge.ls.dat'))

#############################################

# Final spectrum and cluster location plots
pdf(file.path(output_dir, "Cluster_and_spectrum_plots.pdf"), paper="a4", width = 8, height=10)

post.cluster.pos <- post.param.dist_posthoc(gs.out = lymph.gs, ls.merge.out = lymph.merge.ls, 
                                            centiles = c(0.025, 0.5, 0.975), ls.type = "ECR", 
                                            burn.in = number.burn.in, samp.names = short.names, 
                                            plot.heatmap = TRUE, min.threshold = par_min.threshold)  

mut.spec.by.clust <- context.extract_posthoc(mut.spec, lymph.merge.ls, min.threshold = par_min.threshold)

dev.off()

saveRDS(post.cluster.pos, file.path(output_dir, 'object_post.cluster.pos.dat'))

## Write heatmap data to disk
cluster_tbl = tbl_df(post.cluster.pos$heat.dat*2)
cluster_tbl$cluster_id = rownames(post.cluster.pos$heat.dat)

# to handle case of NA clusters in post.param.dist_posthoc
# but NA clusters should disappear by increasing iterations to 16k
# est.num.muts = as.numeric(post.cluster.pos$num.mut[which(post.cluster.pos$num.muts>=par_min.threshold)])
# if(length(post.cluster.pos$idx.na) > 0){
#  est.num.muts = as.numeric(post.cluster.pos$num.mut[which(post.cluster.pos$num.muts>=par_min.threshold)])[-post.cluster.pos$idx.na]
# }

cluster_tbl$estimated.no.of.mutations = as.numeric(post.cluster.pos$num.mut[which(post.cluster.pos$num.muts>=par_min.threshold)])
# cluster_tbl$estimated.no.of.mutations = est.num.muts
cluster_tbl$no.of.mutations.assigned = cluster_tbl$estimated.no.of.mutations

# Add a root node 
heatmap_tbl = rbind(cluster_tbl %>% dplyr::select(-cluster_id),
                    c(rep(1, dim(post.cluster.pos$heat.dat)[2]), sum(cluster_tbl$no.of.mutations.assigned), max(cluster_tbl$no.of.mutations.assigned)))

fwrite(heatmap_tbl, file.path(output_dir, 'heatmap.tsv'), sep='\t')

## Export a table with mutations per cluster
mut_per_cluster = tibble(cluster_id = names(post.cluster.pos$num.mut), num_mut = post.cluster.pos$num.mut)
fwrite(mut_per_cluster, file.path(output_dir, 'muts_per_cluster.csv'))

## Export a heatmap table that includes cluster IDs
fwrite(cluster_tbl %>% dplyr::select(-estimated.no.of.mutations), file.path(output_dir, 'cluster_and_samples.csv'))


#############################################

# Density plots of clusters with at least 100 mutations # edited to 50 muts
pdf(file.path(output_dir, "Raw_mutation_VAF_by_cluster_plots.pdf"), width = 12, height = 15)

which.clusters <- as.double(names(which(table(lymph.merge.ls[["final.cluster"]]) >= par_min.threshold)))

for (i in which.clusters) {
  temp.plot <- cluster.density.plot_posthoc(lymph.gs, lymph.merge.ls, post.cluster.pos, i, samp.names = short.names)
  print(temp.plot)
}

dev.off()

#############################################

# Boxplot of VAFs per sample/cluster using 95.0% centiles (as conf. intervals)
# Added by Fede (fa8):
final_clusters <- unique(ls_tbl$clust_assign)
samples        <- names(post.cluster.pos$centiles[1,,1])
num_muts       <- length(ls_tbl$clust_assign)
colores        <- rainbow(length(final_clusters))

for(sample in samples) {
  pdf(file.path(output_dir, sprintf("Cluster_cellular_prevalence_%s.pdf", sample)), width = 20, height=6)
  #par(mfrow=c(length(samples),1))
  par(mar=c(3,5,1,3))
  plot(NULL,NULL,xlim=c(0,num_muts+length(final_clusters)),ylim=c(0,1.1),xaxt='n',xlab="",ylab="Cell. prev")
  title(sample, line = -1)
  for(h in c(.2,.4,.6,.8)) {
    abline(h=h,col="gray"); 
  }
  xpos  <- 0
  index <- 1
  for(cluster_index in final_clusters) {
    centiles <- post.cluster.pos$centiles[cluster_index,sample,]
    centiles <- centiles * 2 # cellular prevalence
    num_muts_in_cluster <- length(which(ls_tbl$clust_assign==cluster_index))
    #rect(xleft, ybottom, xright, ytop, density = NULL, angle = 45,
    #     col = NA, border = NULL, lty = par("lty"), lwd = par("lwd"),
    #     ...)
    rect(xpos,centiles[1],xpos+num_muts_in_cluster,centiles[3],col=colores[index])
    lines(c(xpos,xpos+num_muts_in_cluster),c(centiles[2],centiles[2]),col="black",lwd=2)
    text((xpos + xpos + num_muts_in_cluster)/2,centiles[2]+0.075,cluster_index)
    xpos  <- xpos  + num_muts_in_cluster + 2
    index <- index + 1 
  }
  dev.off()
}


#save.image()
save.image(file.path(output_dir, 'Rsession_ls2.dat'))

## Write out a file in a format that is convenient to use in downstream analyses:
mut_assignment <- data.frame(Chr=ls_tbl$chrom, Pos=ls_tbl$pos, Ref=lymph.depth$ref,
                             Alt=lymph.depth$alt, Cluster=ls_tbl$clust_assign,Patient=patient,
                             ClusterID=paste(patient, ls_tbl$clust_assign, sep="_"))

write.table(mut_assignment, paste(output_dir, "sbs_assignment.txt", sep=""), sep="\t", row.names = F, quote=F)







