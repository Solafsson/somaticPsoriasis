
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
toutmeans <- read.table("/Users/so11/phd/psoriasis/scratch/branch_exposures_w_prior.txt")
tree_dir="/Users/so11/phd/psoriasis/scratch/trees/"

sample_meta <- read.table("/Users/so11/phd/psoriasis/sample_info/sample_meta.txt", h=T, stringsAsFactors = F)
patients <- unique(sample_meta$patient_ID[sample_meta$Seq_status=="WES" & is.na(sample_meta$excl_criteria)])


## Hard-coded. May need to change if you re-run the HDP
#sigNames <- c("Unassigned","SBS7b - UV","PUVA treatment", "SBS1/5", "UV-component N1",
#              "UV-component N2","SBS2 - APOBEC","SBS7d - UV", "New component N3",
#              "SBS13 - APOBEC")
sigNames <- colnames(toutmeans)
## patients that are still missing: 37,79,
all_sig_muts <- data.frame()
for(patient in patients) {
  
  clusters <- read.csv(paste("/Users/so11/phd/psoriasis/scratch/cluster_and_samples/", patient, "_cluster_and_samples.csv", sep=""), h=T)
  clusters$cluster_id <- gsub("Cl.", paste(patient, "_", sep=""), clusters$cluster_id)
  
  tree <- read.tree(paste(tree_dir, patient, ".cluster_tree.phylo", sep=""))
  tree_original <- tree
  tree_df <- fortify(tree_original)
  rownames(tree_df) <- tree_df$label
  
  ## Extract the branches for this patient. 
  b <- toutmeans[grep(paste(patient, "_", sep=""), rownames(toutmeans)),]
  b$node <- unlist(strsplit(rownames(b), split="_"))[c(FALSE, TRUE)]
  
  test <- merge(tree_df, b, by.x="label", by.y="node", all.x=T)
  # Sort the test dataframe
  rownames(test) <- test$label
  test <- test[rownames(tree_df),]
  
  ## For each signature, multiply the branch lengths with the fraction of mutations attributed to that signature
  ## in that branch. Then sum over all 'anchestor branches' to obtain total counts for each tip. 
  test$label <- paste(patient, test$label, sep="_")
  this_pat_muts <- data.frame(test[test$isTip, c(1,5)])
  for(i in 1:length(sigNames)) {
    
    ## First adjust the length of all branches in the tree by the fraction of mutations.
    adj.branch <- test$branch.length
    adj.branch[!is.na(test[,(9+i)])] <- test$branch.length[!is.na(test[,(9+i)])]*test[!is.na(test[,(9+i)]), (9+i)]
    tree$edge.length=adj.branch[tree$edge[,2]]
    tree$edge.length[is.na(tree$edge.length)] <- 0
    tree_df <- fortify(tree)
    tree_df$label <- paste(patient, tree_df$label, sep="_")
    colnames(tree_df)[4] <- sigNames[i]
    
    ## The exposure in each clone is simply the x-coordinates
    this_pat_muts <- merge(this_pat_muts, tree_df[tree_df$isTip, c(6,4)])
  }
  
  a <- clusters[clusters$cluster_id %in% this_pat_muts$label, grep("P", colnames(clusters))]
  this_pat_muts$HighCellFrac_sample <- colnames(a)[max.col(a,ties.method="first")]
  all_sig_muts <- rbind(all_sig_muts, this_pat_muts)
}  

write.table(all_sig_muts, "/Users/so11/phd/psoriasis/scratch/clone_mutation_burden.txt", sep = "\t", row.names = F, quote = F)
 