
  patientID <- "patient37"
  ndp_input_dir="/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/patient37/ndp_patient37_2021_12_09/"
  script.dir = "/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/06_phylogenics/hdp_clustering/"
  tree_out_dir="/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/trees/"
  min_vaf_threshold=0.05
  min_mut_count=15
  
 .libPaths("/lustre/scratch119/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")
library(RColorBrewer)
library(philentropy)
library(data.table)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(igraph)
library(ape)
library(stringr)
source(paste0(script.dir,"/ndp_tree_functions.R"))
source(paste0(script.dir,"/phylogeny_functions.R"))

repo_location <- paste(script.dir, "snv-clustering-using-the-dirichlet-process/", sep="")

get_cluster_pairs <- function(x.clusters.unique.real.merged, min_num_mut = min_mut_count, min_contrib = min_vaf_threshold){
  target_clust = rownames(x.clusters.unique.real.merged)
  
  heat_cut = x.clusters.unique.real.merged[which(rownames(x.clusters.unique.real.merged) %in% target_clust),]>min_contrib
  
  pair_list = list()
  for(i in 1:dim(heat_cut)[2]){
    inc_clust = names(which(heat_cut[,i]))
      print(inc_clust)
    if(length(inc_clust)<2) { next }
    clust_combs = t(combn(inc_clust, m=2))
    this_pairs = tibble(cl1 = clust_combs[,1], cl2 = clust_combs[,2])
    pair_list[[i]] = this_pairs
  }
  all_pairs = do.call('rbind', pair_list) %>% distinct()
  all_pairs = all_pairs %>%
    mutate(pair_id = sprintf('%s.%s', cl1, cl2), pair_rev_id = sprintf('%s.%s', cl2, cl1))
  all_pairs = all_pairs %>% filter(!(pair_id %in% pair_rev_id))
  
  return(all_pairs %>% dplyr::select(-pair_id, -pair_rev_id))
}


load("/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/patient79/ndp_patient79_2021_12_09/patient79.x.clusters.RData")
load("/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/patient79/ndp_patient79_2021_12_09/patient79.x.clusters.unique.fullrun.RData")
load("/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/patient79/ndp_patient79_2021_12_09/Rsession_ls2.dat")
load("/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/patient79/ndp_patient79_2021_12_09/Rsession_ls.dat")
load("/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/patient79/ndp_patient79_2021_12_09/Rsession.dat")

### plot pigeonhold plots
centiles_tbl_all <- create.diag.plots(x.cluster_and_samples = x.cluster_and_samples.final, x.clusters.centiles = x.clusters.centiles, outdir = ndp_input_dir, id.patient = patientID, num.muts = num.muts, min_num_mut = min_mut_count)

x.cluster_and_samples = x.cluster_and_samples.final
outdir = ndp_input_dir
id.patient = patientID
num.muts = num.muts
min_num_mut = min_mut_count

create.diag.plots
function(x.cluster_and_samples, x.clusters.centiles, outdir, id.patient, num.muts, min_num_mut) {
  target_clust = which(num.muts > min_num_mut)
  out = plot_cluster_vaf(x.cluster_and_samples, x.clusters.centiles, target_clust, outdir, id.patient)
  centiles_tbl_all = na.omit(do.call('rbind', out))
  return(centiles_tbl_all )
}

x.clusters.unique.real.merged=x.cluster_and_samples
x.clusters.unique.real.centiles.final=x.clusters.centiles

> plot_cluster_vaf
function(x.clusters.unique.real.merged, x.clusters.unique.real.centiles.final, target_clust, outdir, id.patient){
  # Generate all cluster pairs
  clust_pairs = get_cluster_pairs(x.clusters.unique.real.merged)
  
  grobs = list()
  centiles_tbl.out = list(data.frame())
  for(i in 1:dim(clust_pairs)[1]){
    this_pair = clust_pairs[i,]
    cl1_id = this_pair$cl1
    cl2_id = this_pair$cl2
    
    message(sprintf('Plotting cluster pair %s and %s', cl1_id, cl2_id))
    
    target_samp = colnames(x.clusters.unique.real.merged)
    
    cluster_centiles = list()
    for(k in c(1:length(target_samp))) {
      idx1 = which(x.clusters.unique.real.centiles.final$id == ifelse(nchar(target_samp[k]) == 15, str_sub(target_samp[k], -6, -1), target_samp[k])   & x.clusters.unique.real.centiles.final$cluster_id == cl1_id)
      idx2 = which(x.clusters.unique.real.centiles.final$id == ifelse(nchar(target_samp[k]) == 15, str_sub(target_samp[k], -6, -1), target_samp[k]) & x.clusters.unique.real.centiles.final$cluster_id == cl2_id)
      if(length(idx1) == 0 || length(idx2) == 0){
        cl1 = rep(NA,3)
        cl2 = rep(NA,3)
        names(cl1) = c('Centile.0.025','Centile.0.5','Centile.0.975')
        names(cl2) = c('Centile.0.025','Centile.0.5','Centile.0.975')
      } else {
        cl1 = x.clusters.unique.real.centiles.final[idx1,c('Centile.0.025','Centile.0.5','Centile.0.975')]
        cl2 = x.clusters.unique.real.centiles.final[idx2,c('Centile.0.025','Centile.0.5','Centile.0.975')]
      }
      
      cluster_centiles[[k]] = data.frame(
        cl1 = this_pair[1],
        cl2 = this_pair[2],
        sample = target_samp[k],
        cl1_Centile.0.025 = as.numeric(cl1['Centile.0.025']),
        cl1_Centile.0.5 = as.numeric(cl1['Centile.0.5']),
        cl1_Centile.0.975 = as.numeric(cl1['Centile.0.975']),
        cl2_Centile.0.025 = as.numeric(cl2['Centile.0.025']),
        cl2_Centile.0.5 = as.numeric(cl2['Centile.0.5']),
        cl2_Centile.0.975 = as.numeric(cl2['Centile.0.975'])
      )
    }
    
    centiles_tbl = tbl_df(as.data.frame((do.call('rbind', cluster_centiles))))
    centiles_tbl.out[[i]] = as.data.frame(centiles_tbl)
    names(centiles_tbl.out)[i] = paste(cl1_id,cl2_id,sep='_')
    
    grobs[[i]] = centiles_tbl %>% 
      ggplot(aes(x=2*cl1_Centile.0.5, y=2*cl2_Centile.0.5)) +
      geom_errorbar(aes(ymin=2*cl2_Centile.0.025, ymax=2*cl2_Centile.0.975), color='red', width=0) +
      geom_errorbarh(aes(xmin=2*cl1_Centile.0.025, xmax=2*cl1_Centile.0.975), color='red', height=0) +
      theme_light() +
      geom_hline(yintercept = 0.5) +
      geom_vline(xintercept = 0.5) +
      geom_abline(intercept = 0, slope=1, size=0.2) +
      geom_abline(intercept = c(0.7,1), slope=-1, size=0.2, color='grey') +
      # geom_text(aes(label = sample), nudge_x = 0.07, nudge_y = 0.025, size=1) +
      geom_text(aes(label = sample), position=position_jitter(width=0.025,height=0.025), size=1) +
      scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
      scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
      labs(x=sprintf('Cell fraction, cluster %s.', cl1_id), y=sprintf('Cell fraction, cluster %s.', cl2_id))
     ggsave(file.path(outdir, sprintf('pigeonhole_clusters_%s_%s.pdf', cl1_id, cl2_id)), width=12, height=12, units='cm')
  }
  # Plot all cluster pairs in a single figure
  pdf(file.path(outdir, sprintf('pigeonhole_clusterpairs_%s.pdf', id.patient)), w=20, h=20)
  do.call('grid.arrange', grobs)
  dev.off()
  return(centiles_tbl.out)
}



muts_per_cluster <- tibble(mut_df)
x.branches.keep.final <- tibble(From="P", To=mut_df$cluster_id, Evidence="strong")
tmp <- mut_df
rownames(tmp) <- mut_df$cluster_id

### create a new ape tree object
tr = convert_edges_to_phylo(x.branches.keep.final, muts_per_cluster)
tr$edge.length <- tmp[tr$tip.label, "num.muts"]

### save treefile
write.tree(tr, file = paste(tree_out_dir,paste(patientID, 'cluster_tree.phylo',sep='.'),sep='/'))

### plot tree
pdf(paste(tree_out_dir,paste(patientID, 'cluster_tree.pdf',sep='.'),sep='/'))
par(xpd=TRUE, mar=c(12.1, 4.1, 4.1, 2.1), family = "sans")
draw_nice_tree(tr,show_internal_nodes=T,edge_strength=x.branches.keep.final$Evidence, patientID=patientID)
highlight_tree_nodes(tr, pch=23, bg='red',show_internal_nodes=T)
dev.off()
