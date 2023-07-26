# Script for generation trees from NDP output

#############################################
## -- Receive command line args --
args = commandArgs(trailingOnly=TRUE)
args

if(FALSE) {
  patientID <- "patient01"
  ndp_input_dir=paste("/lustre/scratch126/humgen/projects/psoriasis/phylogenics/hdp_clustering/", patientID, "/ndp_", patientID, "_2021_12_09/", sep="")
  script.dir = "/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/06_phylogenics/hdp_clustering/"
  tree_out_dir="/lustre/scratch126/humgen/projects/psoriasis/phylogenics/hdp_clustering/trees_updated/"
  min_vaf_threshold=0.05
  min_mut_count=10
}

# Check for command line args, else set defaults
if (length(args)>=3) {
  patientID = args[1]
  ndp_input_dir = args[2]
  script.dir = args[3]
} else {
  stop("Script needs at least 4 arguments.", call.=FALSE)
} 

if (length(args) >= 4) {
  tree_out_dir = args[4]
} else { 
  tree_out_dir=ndp_input_dir
}	

if (length(args) >= 5) {
  min_vaf_threshold = args[5]
} else { 
  min_vaf_threshold =0.10
  message("Using default min VAF threshold of 0.10")
}	

if (length(args) >= 6) {
  min_mut_count = args[6]
} else { 
  min_mut_count =50
  message("Using default min mutation threshold of 50")
}	

snv_dir =ndp_input_dir
#cat(paste(snv_dir,paste0(patientID,'_bb_pass_snvs_all.csv'),sep='/'))

### Load required Libraries

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
source(paste0(script.dir,"/ndp_tree_functions.R"))
source(paste0(script.dir,"/phylogeny_functions.R"))


repo_location <- paste(script.dir, "snv-clustering-using-the-dirichlet-process/", sep="")

### Get clean mut spectrum
mut.spec.by.clust <- extract.mut.spec(ndp_input_dir)

### Assign mut cluser ID
x.cluster_and_samples <- cluster.append(snvdir =  snv_dir, ndpdir = ndp_input_dir,  mut.spec.by.clust = mut.spec.by.clust, id.patient = patientID)

### Filter to remove clusters not present in any sample above min vaf threshold
x.cluster_and_samples <- x.cluster_and_samples %>% filter(apply(x.cluster_and_samples[, c(1:(ncol(..x.cluster_and_samples)-2))], 1, function(x) any (x >= min_vaf_threshold)))

#save(x.cluster_and_samples, mut.spec.by.clust,file=paste(ndp_input_dir,paste(patientID,'x.clusters.RData',sep='.'),sep='/'))

### get number of mutations for branch lengths
mut_df <- extract.mutation.counts(x.cluster_and_samples = x.cluster_and_samples)
num.muts = x.cluster_and_samples$no.of.mutations.assigned
names(num.muts) = x.cluster_and_samples$cluster_id

### extract 95% confidence interval data
x.clusters.centiles <- extract.95.CI(x.cluster_and_samples = x.cluster_and_samples, ndpdir = ndp_input_dir)

### format median VAF tbl
x.cluster_and_samples.final <- create.med.vaf.tbl(x.cluster_and_samples = x.cluster_and_samples)

### asssign clusters to annotated mutations
file.snv = paste(snv_dir,paste0(patientID,'_bb_pass_snvs_all.csv'),sep='/')
snv_tbl <- fread(file.snv) %>%
  dplyr::mutate(pos_id = paste(chrom, pos, sep = "_"))
cluster.file <- paste(ndp_input_dir, "clust_assign_posthoc.csv", sep = "/")
cluster_tbl <- fread(cluster.file) %>%
  dplyr::mutate(pos_id = paste(chrom, pos, sep = "_")) %>%
  dplyr::select(-mut_id, -chrom, -pos)
assigned_tbl <- left_join(snv_tbl, cluster_tbl, by = "pos_id") %>%
  dplyr::rename("cluster_id" = clust_assign)
#write.csv(assigned_tbl, file = paste(snv_dir, paste0(patientID, "_ndp_assigned_muts.csv"), sep = "/"), quote = F, row.names = F)

### save intermediate r data
#save(mut_df, num.muts, x.clusters.centiles,x.cluster_and_samples.final,file=paste(ndp_input_dir,paste(patientID,'x.clusters.unique.fullrun.RData',sep='.'),sep='/'))

### generate all combinations of cluster pairs that pass thresholds
get_cluster_pairs <- function(x.clusters.unique.real.merged, min_num_mut = min_mut_count, min_contrib = min_vaf_threshold){
  target_clust = rownames(x.clusters.unique.real.merged)
  
  heat_cut = x.clusters.unique.real.merged[which(rownames(x.clusters.unique.real.merged) %in% target_clust),]>min_contrib
  
  pair_list = list()
  for(i in 1:dim(heat_cut)[2]){
    inc_clust = names(which(heat_cut[,i]))
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

### plot pigeonhold plots
centiles_tbl_all <- create.diag.plots(x.cluster_and_samples = x.cluster_and_samples.final, x.clusters.centiles = x.clusters.centiles, 
                  outdir = ndp_input_dir, id.patient = patientID, num.muts = num.muts, min_num_mut = min_mut_count) 
clust_pairs = unique(centiles_tbl_all[,c('cl1','cl2')])

### determine tree branches 
x.branches <-  infer.tree.branches(clust_pairs = clust_pairs, centiles_tbl_all = centiles_tbl_all, min_contrib = min_vaf_threshold)

## Add "Independent cluster branches"
## Added by Sigurgeir
idx = which(rownames(x.cluster_and_samples.final) %in% rownames(x.cluster_and_samples.final)[!(rownames(x.cluster_and_samples.final) %in% x.branches)])
if(length(idx)>0) {
  for(i in 1:length(idx)){
    x.branches = rbind(x.branches,c('P',rownames(x.cluster_and_samples.final)[idx[i]],'strong'))
  }
}



# decide which branches to keep
x.branches.keep = keep(x.branches)
x.branches.keep = x.branches.keep[!duplicated(x.branches.keep[,c('From','To')]),]
idx = which(x.branches.keep$From %in% 'P')
ids.clusters = unique(sort(as.character(x.branches.keep$From)[-idx]))
idx = which(ids.clusters %in% x.branches.keep$To[-idx])
if(length(idx) > 0){
  ids.clusters = ids.clusters[-idx]
}
x.branches.keep$Evidence = NA
rownames(x.branches.keep) = seq(1,dim(x.branches.keep)[1])

# assign strength of evidence of clonal relationships
x.branches = as.data.frame(x.branches)
for(i in 1:dim(x.branches.keep)[1]){
  idx = which(x.branches$From %in% x.branches.keep$From[i] & x.branches$To %in% x.branches.keep$To[i])
  x.branches.keep$Evidence[i] = as.character(x.branches$Evidence[idx])
}
#save(x.branches.keep, centiles_tbl_all,file=paste(ndp_input_dir,paste(patientID,'x.branches.pre_tiebreak.RData',sep='.'),sep='/'))

### break ties if multiple potential paths, keeping higher VAF
x.branches.keep.final <- break.ties(x.branches.keep = x.branches.keep, centiles_tbl_all = centiles_tbl_all)

### save intermediate r data
#save(mut_df, x.branches.keep.final,file=paste(ndp_input_dir,paste(patientID,'x.branches.keep.RData',sep='.'),sep='/'))

### prep tree building data
x.branches.keep.final$To = gsub('Cl.','',x.branches.keep.final$To,ignore.case=T)
x.branches.keep.final$From = gsub('Cl.','',x.branches.keep.final$From,ignore.case=T)
muts_per_cluster = mut_df

### create a new ape tree object
tr = convert_edges_to_phylo(x.branches.keep.final, muts_per_cluster)


### save treefile
write.tree(tr, file = paste(tree_out_dir,paste(patientID, 'cluster_tree.phylo',sep='.'),sep='/'))

### plot tree
pdf(paste(tree_out_dir,paste(patientID, 'cluster_tree_newEvidence.pdf',sep='.'),sep='/'))
par(xpd=TRUE, mar=c(12.1, 4.1, 4.1, 2.1), family = "sans")
draw_nice_tree(tr,show_internal_nodes=T,edge_strength=x.branches.keep.final$Evidence[c(nrow(x.branches.keep.final):1)], patientID=patientID)
highlight_tree_nodes(tr, pch=23, bg='red',show_internal_nodes=T)
dev.off()
