
## The purpose of this tree is to prune from the phylogenetic trees all branches where the 
## pigeonhole principle may not be valid. We want to see if our estimate of the mutation rate
## varies if we include only clusters where the pigeonhole is incontrovertible. 


.libPaths("/lustre/scratch126/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")
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
script.dir = "/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/06_phylogenics/hdp_clustering/"
source(paste0(script.dir,"/ndp_tree_functions.R"))
source(paste0(script.dir,"/phylogeny_functions.R"))

patient_meta <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/Supplementary_material/Supplementary_Table1_patient_metadata.txt", h=T, sep = "\t")

## For each patient, I want to do the following: 
## 1. Write out a plot of the full tree with the correct branch evidence
## 2. Prune the tree, retaining only tips with "strong" evidence from the root to the tip or tips where a strong branch follows a weaker branch. 
## 3. Write out the tree and plot it. 
## 4. (in a different script) estimate the mutational signatures and mutation burden of each tree tip (clone) as before

for(patientID in patient_meta$Patient.ID) {
  if(patientID %in% c("patient79")) {
    next
    ndp_input_dir=paste("/lustre/scratch126/humgen/projects/psoriasis/phylogenics/hdp_clustering/", patientID, "/ndp_", patientID, "_2022_02_15/", sep="")
    
  } else {
    ndp_input_dir=paste("/lustre/scratch126/humgen/projects/psoriasis/phylogenics/hdp_clustering/", patientID, "/ndp_", patientID, "_2021_12_09/", sep="")
    
  }
  tree_out_dir=ndp_input_dir
  snv_dir =ndp_input_dir
  repo_location <- paste(script.dir, "snv-clustering-using-the-dirichlet-process/", sep="")
  new_tree_dir="/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/pruned_trees/"
  
  ### Get clean mut spectrum
  #mut.spec.by.clust <- extract.mut.spec(ndp_input_dir)
  ### Assign mut cluser ID
  #x.cluster_and_samples <- cluster.append(snvdir =  snv_dir, ndpdir = ndp_input_dir,  mut.spec.by.clust = mut.spec.by.clust, id.patient = patientID)
  
  load(paste(ndp_input_dir,paste(patientID,'x.branches.keep.RData',sep='.'),sep='/'))
  
  ### prep tree building data
  x.branches.keep.final$To = gsub('Cl.','',x.branches.keep.final$To,ignore.case=T)
  x.branches.keep.final$From = gsub('Cl.','',x.branches.keep.final$From,ignore.case=T)
  muts_per_cluster = mut_df
  
  if(patientID=="patient37") {
    x.branches.keep.final <- x.branches.keep.final[c(1:2,4:10),]
  }
  
  ### create a new ape tree object
  tr = convert_edges_to_phylo(x.branches.keep.final, muts_per_cluster)
  
  ## Decide which branches to keep
 # x.branches.keep.final$keep <- F
 # x.branches.keep.final$keep[x.branches.keep.final$Evidence=="strong"] <- T
 # x.branches.keep.final$keep[]
  
  
  x.branches.pruned <- x.branches.keep.final[x.branches.keep.final$Evidence=="strong",]
  x.branches.pruned$ID <- paste("Cl.", x.branches.pruned$To, sep="")
  muts_pruned <- muts_per_cluster[muts_per_cluster$cluster_id %in% x.branches.pruned$ID,]
  pruned_tree <- convert_edges_to_phylo(x.branches.pruned, muts_pruned)
  
  ### save treefile
  write.tree(tr, file = paste(new_tree_dir,paste(patientID, 'cluster_tree.phylo',sep='.'),sep='/'))
  write.tree(pruned_tree, file = paste(new_tree_dir,paste(patientID, 'pruned_cluster_tree.phylo',sep='.'),sep='/'))

  ### plot tree
  pdf(paste(new_tree_dir,paste(patientID, 'cluster_tree.pdf',sep='.'),sep='/'))
  par(xpd=TRUE, mar=c(12.1, 4.1, 4.1, 2.1), family = "sans")
  draw_nice_tree(tr,show_internal_nodes=T,edge_strength=x.branches.keep.final$Evidence[c(nrow(x.branches.keep.final):1)], patientID=patientID)
  #highlight_tree_nodes(tr, pch=23, bg='red',show_internal_nodes=T)
  dev.off()
  
  ### plot tree
  pdf(paste(new_tree_dir,paste(patientID, 'pruned_cluster_tree.pdf',sep='.'),sep='/'))
  par(xpd=TRUE, mar=c(12.1, 4.1, 4.1, 2.1), family = "sans")
  draw_nice_tree(pruned_tree,show_internal_nodes=T,edge_strength=x.branches.pruned$Evidence[c(nrow(x.branches.pruned):1)], patientID=patientID)
  #highlight_tree_nodes(tr, pch=23, bg='red',show_internal_nodes=T)
  dev.off()
  
}

