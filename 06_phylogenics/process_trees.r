## This script does the following:
## 1. Read in a MPBoot tree, genotype matrix and the nr of reads supporting mutation/ref alleles.
## 2. Use maximum Likelihood to assing each mutation to a branch of the tree
## 3. Plot the tree after this
## 4. Adjust the branch lengths by VAF and coverage. Plot a new tree and sum the mutation counts for
## each sample to be used in lmm of mutation burden.
## 5. Assign mutations to branches and write this out for signature extraction and driver analyses. 

#.libPaths("/lustre/scratch114/projects/crohns/somatic_ibd_p1/R_packages_farm5_R3.6_install/")
.libPaths("/lustre/scratch119/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")
library(ape)
library(ggtree)
library(VGAM)
library(ggplot2)

source("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/06_phylogenics/treemut.R")


## Define parameters
if(T) {
  
  sample_meta <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/sample_info/sample_meta.txt", h=T)
  
  parameters <- commandArgs(TRUE)
  PID <- parameters[1]
  gt_matrix <- parameters[2]
  treeFile <- parameters[3]
  NR_file <- parameters[4]
  NV_file <- parameters[5]
  output_dir <- parameters[6]
  mutType <- parameters[7]
  

  
  
} else {
  
  ## For locally running and testing: 
  sample_meta <- read.table("/Users/so11/phd/psoriasis/sample_info/sample_meta_wrongNames_first148.txt", h=T)
  PID <- "patient01"
  gt_matrix <- "/Users/so11/phd/psoriasis/scratch/binomial_filters/patient01/patient01_genotype_sbs.txt"
  treeFile <- "/Users/so11/phd/psoriasis/scratch/phylogenics/MPBoot/snv/consensus_trees/patient01/patient01.treefile"
  NR_file <- "/Users/so11/phd/psoriasis/scratch/binomial_filters/patient01/patient01_NR_pass_sbs.txt"
  NV_file <- "/Users/so11/phd/psoriasis/scratch/binomial_filters/patient01/patient01_NV_pass_sbs.txt"
  output_dir <- "/Users/so11/phd/psoriasis/scratch/"
  mutType <- "snv"
    
  source("/Users/so11/phd/psoriasis/scratch/treemut.R")
  library(ape)
  library(ggtree)
  library(VGAM)
}

options(stringsAsFactors = F)

## Read in files: 
NR <- read.table(NR_file, h=T)
NV <- read.table(NV_file, h=T)
Genotype <- read.table(gt_matrix, h=T)
tree=read.tree(treeFile)


## Clean the tree as well. 
tree$tip.label=gsub("_VAF","",tree$tip.label)
tree$edge.length=rep(1,nrow(tree$edge))
tree=drop.tip(tree,"Ancestral")
NR_flt = as.matrix(NR[rownames(Genotype),tree$tip.label])
NV_flt = as.matrix(NV[rownames(Genotype),tree$tip.label])

## Get the meta-data sorted:
meta <- sample_meta[, c("sampleID", "biopsy_ID", "Location", "type")]



## Assign mutations to the tree
df=reconstruct_genotype_summary(tree)
res=assign_to_tree(df=df,
                   mtr=NV_flt,
                   dep=NR_flt)


## Make the edge lengths equal 
edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<0.01])
#edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<=1])
edge_length = rep(0,nrow(tree$edge))
names(edge_length)=1:nrow(tree$edge)
edge_length[names(edge_length_nonzero)]=edge_length_nonzero

tree2=tree
tree2$edge.length=as.numeric(edge_length)

system(paste("mkdir -p ", output_dir, "consensus_trees/", PID, "/", sep=""))
system(paste("mkdir -p ", output_dir, "tree_plots/", PID, "/", sep=""))

write.tree(tree2, paste(output_dir, "consensus_trees/", PID, "/", PID, "_MutAssign_unadj.tree", sep=""))


pdf(paste(output_dir, "tree_plots/", PID, "/", PID, "_MutAssign_unadj.pdf", sep=""), onefile=T)
plot <- ggtree(tree2) 
plot <- plot %<+% meta + geom_tiplab(aes(color=type),vjust=-0.3, size=6) +
  geom_tippoint(aes(color=type), alpha=0.25) + 
scale_color_manual(values=c(Healthy = "#FFDBAC",Lesional= "#e18377"))
plot <- plot +theme_tree2()+xlim(0,max(fortify(tree2)$x)*1.3) + 
  theme(legend.position="right", legend.title = element_blank(), legend.text = element_text(size=16))
print(plot)


plot <- ggtree(tree2) 
plot <- plot %<+% meta + geom_tiplab(aes(color=type),vjust=-0.3, size=6) +
  geom_tippoint(aes(color=type), alpha=0.25) + 
scale_color_manual(values=c(Healthy = "#FFDBAC",Lesional= "#e18377"))
plot <- plot +theme_tree2()+xlim(0,max(fortify(tree2)$x)*1.3) + 
  geom_text2(aes(label=node), hjust=-.3, vjust=0.9, color="red") + 
  theme(legend.position="right", legend.title = element_blank(), legend.text = element_text(size=16))
print(plot)
dev.off()


## Write out a file assigning mutations to branches. 
system(paste("mkdir -p ", output_dir, "branch_mut_assignment/", PID, "/", sep=""))

con = textConnection(rownames(Genotype))
Mutations_per_branch = read.table(con,sep=":")
colnames(Mutations_per_branch)=c("Chr","Pos","Ref","Alt")
Mutations_per_branch$Branch = tree$edge[res$summary$edge_ml,2]
Mutations_per_branch$Branch[rowSums(Genotype>0.3)==ncol(Genotype)]=0 #Root
Mutations_per_branch=Mutations_per_branch[Mutations_per_branch$Branch!=0&res$summary$p_else_where<0.01,]
Mutations_per_branch$Patient = PID
Mutations_per_branch$SampleID = paste(PID,Mutations_per_branch$Branch,sep="_")
write.table(Mutations_per_branch,paste(output_dir, "branch_mut_assignment/", PID, "/",PID, "_ML_branch_assignment.txt", sep=""),quote=F,row.names=F,sep="\t")



## Need to figure out if I'm going to adjust the mutation count for coverage and VAF of the sample
## See functions in treemut.R
