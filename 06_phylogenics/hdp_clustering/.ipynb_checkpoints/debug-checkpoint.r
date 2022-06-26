
  patientID <- "patient79"
  ndp_input_dir="/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/patient79/ndp_patient79_2021_12_09/"
  script.dir = "/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/06_phylogenics/hdp_clustering/"
  tree_out_dir="/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/trees/"
  min_vaf_threshold=0.05
  min_mut_count=10
  
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


