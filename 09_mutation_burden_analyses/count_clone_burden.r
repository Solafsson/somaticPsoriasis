
## The purpose of this script is to create new "pseudo-samples" or clones for each tip of the phylogenetic
## tree of each patient. I also want to estimate the contribution of each signature to the total mutation
## burden of each clone. 

#bsub -o count_mutations.out -e count_mutations.err -R 'select[mem>=1000] rusage[mem=1000]' -M1000  "/software/R-4.1.0/bin/Rscript /nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/09_mutation_burden_analyses/count_clone_burden.r" 

.libPaths("/lustre/scratch119/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")
library(ape)
library(ggtree)
library(reshape2)
library(ggplot2)
library(cowplot)


#tree_dir="/Users/so11/phd/psoriasis/scratch/phylogenics/MPBoot/snv/consensus_trees/"
#toutmeans <- read.table("/Users/so11/phd/psoriasis/scratch/snv/branch_exposures_w_prior.txt", h=T)
# toutmeans <- read.table("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/hdp/snv/branch_exposures_w_prior.txt", h=T)
toutmeans <- read.table("/Users/so11/phd/psoriasis/scratch/branch_exposures_w_prior.txt")
tree_dir="/Users/so11/phd/psoriasis/scratch/trees/"

sample_meta <- read.table("/Users/so11/phd/psoriasis/sample_info/sample_meta.txt", h=T, stringsAsFactors = F)
patients <- unique(sample_meta$patient_ID[sample_meta$Seq_status=="WES" & is.na(sample_meta$excl_criteria)])


## Hard-coded. May need to change if you re-run the HDP
#sigNames <- c("Unassigned","SBS7b - UV","PUVA treatment", "SBS1/5", "UV-component N1",
#              "UV-component N2","SBS2 - APOBEC","SBS7d - UV", "New component N3",
#              "SBS13 - APOBEC")
sigNames = colnames(toutmeans) 
all_sig_muts <- data.frame()
for(patient in patients) {
  
  clusters <- read.csv(paste("/Users/so11/phd/psoriasis/scratch/cluster_and_samples/", patient, "_cluster_and_samples.csv", sep=""), h=T)
  clusters$cluster_id <- gsub("Cl.", paste(patient, "_", sep=""), clusters$cluster_id)
  
  tree <- read.tree(paste(tree_dir, patient, ".cluster_tree.phylo", sep=""))
  tree$tip.label <- gsub("Cl.", "", tree$tip.label)
  tree_df <- fortify(tree)
  tree_df$label <- paste(patient, tree_df$label, sep="_")
  
  ## read in the corrected branch lengths:
  correction <- read.table(paste("/Users/so11/phd/psoriasis/scratch/cluster_sensitivity/", patient, "_cluster_sensitivity.txt", sep=""), h=T)
  test <- merge(tree_df, correction[,c("ClusterID", "Mutations_adj")], by.x="label", by.y="ClusterID", all.x=T)
  test$branch.length <- test$Mutations_adj
  test$Mutations_adj <- NULL
  test <- test[order(test$node),]
  tree$edge.length <- test$branch.length[tree$edge[,2]]
  
  tree_df <- fortify(tree)
  tree_df$label <- paste(patient, tree_df$label, sep="_")
  ## Extract the clones for this patient. 
  b <- toutmeans[grep(paste(patient, "_", sep=""), rownames(toutmeans)),]
  b$node <- rownames(b)
  
  test <- merge(tree_df, b, by.x="label", by.y="node", all.x=T)
  # Sort the test dataframe
  test <- test[order(test$node),]
  
  ## For each signature, multiply the branch lengths with the fraction of mutations attributed to that signature
  ## in that branch. Then sum over all 'anchestor branches' to obtain total counts for each tip. 
  ## Clusters that have fewer than 20 mutations were not included in the signature extraction.
  ## assign all mutations to the "Unassigned" component and set the rest to zero.
  this_pat_muts <- data.frame(test[test$isTip, c(1,5)])
  test[is.na(test[,10]),c(10:ncol(test))] <- 0
  
  for(i in 1:length(sigNames)) {
    
    ## First adjust the length of all branches in the tree by the fraction of mutations.
    adj.branch <- test$branch.length
    adj.branch[!is.na(test[,(9+i)])] <- test$branch.length*test[, (9+i)]
    tree$edge.length=adj.branch[tree$edge[,2]]
    tree$edge.length[is.na(tree$edge.length)] <- 0
    tree_df <- fortify(tree)
    tree_df$label <- paste(patient, tree_df$label, sep="_")
    colnames(tree_df)[4] <- sigNames[i]
    
    ## The exposure in each clone is simply the x-coordinates
    this_pat_muts <- merge(this_pat_muts, tree_df[tree_df$isTip, c(6,4)])
  }
  
  all_sig_muts <- rbind(all_sig_muts, this_pat_muts)
} 

## Find the sample which in which the clone has the highest cellular prevalence. 
## This will help link the clones with the rest of the meta-data. 

clone_burden <- all_sig_muts
clone_burden$HighCellFrac_sample <- NA
sample_meta <- read.table("/Users/so11/phd/psoriasis/sample_info/sample_meta.txt", h=T, stringsAsFactors = F)
patients <- unique(sample_meta$patient_ID[sample_meta$Seq_status=="WES" & is.na(sample_meta$excl_criteria)])
for(patient in patients) {
  clusters <- read.csv(paste("/Users/so11/phd/psoriasis/scratch/cluster_and_samples/", patient, "_cluster_and_samples.csv", sep=""), h=T)
  clusters$cluster_id <- gsub("Cl.", paste(patient, "_", sep=""),clusters$cluster_id)
  clusters$maxExp <- apply(clusters[,c(1:(ncol(clusters)-2))], 1, max)
  
  for(i in 1:nrow(clusters)) {
    if(clusters$cluster_id[i] %in% clone_burden$label) {
      
      ## This is just helpful for annotation. Which clones come from lesional vs non-lesional biopsies etc. 
      clone_burden$HighCellFrac_sample[clone_burden$label==clusters$cluster_id[i]] <- names(which.max(clusters[i, c(1:(ncol(clusters)-3))]))
      
    }
  }
  
}

colnames(clone_burden)[1:2] <- c("CloneID", "TotalSBS_adj")
write.table(clone_burden, "/Users/so11/phd/psoriasis/scratch/clone_mutation_burden.txt", sep="\t", row.names = F, quote=F)

