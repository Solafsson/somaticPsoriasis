
## multiply the cellular fraction of the cluster and the volume of the sample. 
## Divide by the total volume of the samples from that patient. 
## Assume that all the mutations assigned to a cluster are present in a sample 
## at the cellular prevalence of a cluster. This may mean that we assign a mutation
## to a sample in which it is not called, if other mutations in the cluster are 
## called in the sample. 

## Assign indels and DBS mutations to clusters based on their VAF vectors and the
## similarity of the VAF vectors to the cluster prevalence. I only do this for drivers.

## Start by extracting the cluster for all snvs. 
## for all non-SNVs, find the cluster where the euclidian distance between 
## half of the estimated cell fractions and the VAF-vector for the mutation is
## minimized. Assign the indel/DBS to that cluster. 

#######################################

euclidean <- function(a, b) sqrt(sum((a - b)^2))
mutation_list <- read.table("/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/annotated_coding_mutations.txt", h=T)
mutation_list$mutID <- paste("chr", mutation_list$chr, ":", mutation_list$pos, ":", mutation_list$ref, ":", mutation_list$mut, sep="")
mutation_list$patient <- unlist(strsplit(substr(mutation_list$sampleID, 1,6), split="L|H"))[c(T,F)]
mutation_list$patient <-gsub("P", "patient", mutation_list$patient)

## Extract potential driver mutations that are within significant genes. 
#geneList <- read.table("/lustre/scratch119/humgen/projects/psoriasis/selectionQTLs/genelist.txt")
geneList=data.frame(V1=c("EEF1A1","NOTHC1", "NOTCH2", "TP53", "PPM1D", "FAT1", "CHEK2", "GXYLT1", "ZFP36L2"))

muts <- mutation_list[mutation_list$gene %in% geneList$V1 & mutation_list$impact!="Synonymous", ]

patient_list <- unique(muts$patient)

## Assign all potential driver mutations to clusters. 
## For non-snv mutations, assign the mutations to clusters based on Euklidian distance
## between the cluster cell fractions and the VAF vectors for that mutation. 
mut_clust_assigned <- data.frame()
for(patient in patient_list) {
  
  thisPatient <- muts[muts$patient==patient,]
  if(dim(thisPatient)[1]==0) {
    ## No driver mutations found for this patient
    next
  }
  
  VAFs <- read.table(paste("/lustre/scratch119/humgen/projects/psoriasis/binomial_filters/", patient, "/", patient, "_VAFs_pass_allMuts.txt", sep=""), h=T)
  VAFs$mutID <- rownames(VAFs)
  if(patient=="patient02") {
    ## Hacky - No drivers identified in these two samples anyway...
    VAFs$P02L_11 <- VAFs$P02L_12 <- NULL
  }
  
  cluster_and_samples <- read.csv(paste("/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/", patient, "/ndp_",patient, "_2021_12_09/cluster_and_samples.csv",sep=""), h=T)
  cluster_and_samples$cluster_id <- gsub("Cl.", paste(patient, "_", sep=""), cluster_and_samples$cluster_id)
  sbs_assign <- read.table(paste("/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/", patient, "/ndp_", patient, "_2021_12_09sbs_assignment.txt", sep=""),h=T)
  sbs_assign$mutID <- paste(sbs_assign$Chr, sbs_assign$Pos, sbs_assign$Ref, sbs_assign$Alt, sep=":")
  
  ## For snvs, simply extract the clusterID from the sbs_assign file
  thisPatient <- merge(thisPatient, sbs_assign[,c("mutID", "ClusterID")], by="mutID", all.x=T)
  
  ## for the other mutations, assign them the cluster that has the lowest Euclidian
  ## distance to the VAF vector of the mutation. 
  for(mutation in thisPatient$mutID[is.na(thisPatient$ClusterID)]) {
    
    tmp <- as.numeric()
    for(i in 1:(nrow(cluster_and_samples))) {
      tmp[i] <- euclidean(cluster_and_samples[i,c(1:(ncol(cluster_and_samples)-2))],
                           VAFs[mutation,c(1:(ncol(VAFs)-1))]*2)
    }
    cl <- cluster_and_samples[which(tmp==min(tmp)),c("cluster_id")]
    thisPatient$ClusterID[thisPatient$mutID==mutation] <- cl
  }
  mut_clust_assigned <- rbind(mut_clust_assigned, thisPatient)
}

write.table(mut_clust_assigned, file="/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/mut_clust_assigned.txt", sep="\t", row.names=F, quote=F)


## Next task is to calculate the fraction of cells that carry a mutation in each gene.
## To do this, I will weigh the cellular prevalences by the volume of the sample but
## also do the calculation in a phylogenetically-informed way. If two clusters have
## a mutation in the same gene, but one cluster is a sub-clone of the other, then I want to
## count only the parent clone (which has larger VAF)
options(stringsAsFactors = F)
mut_clust_assigned <- read.table("/Users/so11/phd/psoriasis/scratch/mut_clust_assigned.txt", h=T)
library(ape)
library(ggtree)
cell_frac_input <- mut_clust_assigned
cell_frac_input$isSubclone <- F
cell_frac_input$uniqueID <- paste(cell_frac_input$patient, cell_frac_input$mutID, sep = "_")
patient_list <- as.character(unique(cell_frac_input$patient))
geneList <- data.frame(as.character(unique(cell_frac_input$gene)))
colnames(geneList) <- "V1"

#tree_dir="/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/trees/"
tree_dir="/Users/so11/phd/psoriasis/scratch/trees/"
for(patient in patient_list) {
  thisPat <- cell_frac_input[cell_frac_input$patient==patient,]
  tree <- read.tree(paste(tree_dir, patient, ".cluster_tree.phylo", sep=""))
  for(gene in geneList$V1) {
    thisGene <- thisPat[thisPat$gene==gene,]
    if(dim(thisGene)[1] < 2) {
      # 0 or 1 mutation in the gene for this patient, no need to do anything.
      next
    } else {
      ## Going to mark mutations as belonging to a sub-clone if they have a parent clone in the tree
      ## that has a mutation in the same gene. The cell fraction of the parent clone will be used in the calculation. 
      tree_df <- fortify(tree)
      tree_df$label <- paste(patient, tree_df$label, sep="_")
      for(i in 1:nrow(thisGene)) {
        thisMutation <- thisGene$uniqueID[i]
        cluster <- thisGene$ClusterID[i]
        if(!(cluster %in% tree_df$label)) {
          ## The cluster is not in the tree for some reason
          print(paste(cluster, "is not in the tree"))
          next
        }
        parent <- tree_df$label[tree_df$parent[tree_df$label==cluster]]
        ## patient_P is the root of the tree
        while(parent!= paste(patient, "P", sep="_")) {
          ## Walk down the tree and see if any of the parent clusters are in the mutation list. If so,
          ## set cell_frac_input$isSubclone as True.
          cluster <- parent
          if(cluster %in% thisGene$ClusterID) {
            cell_frac_input$isSubclone[cell_frac_input$uniqueID==thisMutation] <- T
          }
          parent <- tree_df$label[tree_df$parent[tree_df$label==cluster]]
        }
      }
    }
  }
} 

write.table(cell_frac_input, file="/Users/so11/phd/psoriasis/selection_analyses/probable_drivers.txt", sep="\t", row.names = F, quote=F)


##############################################
##
##  CALCULATE FRACTION OF MUTATED CELLS
##############################################

## The final step is to weigh the samples by their volume
## and sum up the cellular prevalence 
geneList <- read.table("/lustre/scratch119/humgen/projects/psoriasis/selectionQTLs/genelist.txt")
cell_frac_input <- read.table("/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/probable_drivers.txt", h=T)
microb_meta <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/Microdissection_metaData.txt", h=T)
patient_list <- as.character(unique(microb_meta$PatientID[microb_meta$ExclusionCriteria=="PASS"]))

cell_frac_results_total <- matrix(nrow=length(patient_list), ncol=length(geneList$V1))
cell_frac_results_lesional <- matrix(nrow=length(patient_list), ncol=length(geneList$V1))
cell_frac_results_non_lesional <- matrix(nrow=length(patient_list), ncol=length(geneList$V1))
r <- 1
for(patient in patient_list) {
  cluster_and_samples <- read.csv(paste("/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/", patient, "/ndp_",patient, "_2021_12_09/cluster_and_samples.csv",sep=""), h=T)
  cluster_and_samples$cluster_id <- gsub("Cl.", paste(patient, "_", sep=""), cluster_and_samples$cluster_id)
  volumes <- microb_meta[microb_meta$PatientID==patient & microb_meta$ExclusionCriteria %in% c("PASS","Technical_duplicate"), c("SampleID", "CutVolume")]
  clusterVolumes <- matrix(nrow = nrow(cluster_and_samples), ncol = nrow(volumes))
  j <- 1
  for(sample in colnames(cluster_and_samples)[grep("P", colnames(cluster_and_samples))]) {
    for(i in 1:nrow(cluster_and_samples)) {
      ## The volume of the LCM cut times the cellular prevalence of the cluster
      ## in each microbiopsy.
      clusterVolumes[i,j] <- cluster_and_samples[i, sample]*volumes$CutVolume[volumes$SampleID==sample]
    }
    j <- j+1
  }
  cluster_and_samples$Volume <- rowSums(clusterVolumes, na.rm=T)
  cluster_and_samples$Volume_lesional <- rowSums(clusterVolumes[,grep("H", colnames(cluster_and_samples)[c(1:(ncol(cluster_and_samples)-4))], invert = T)])
  
  nrNonLesional <- length(grep("H", colnames(cluster_and_samples)))
  if(nrNonLesional==0) {
    cluster_and_samples$Volume_non_lesional <- NA
  } else if(nrNonLesional==1) {
    cluster_and_samples$Volume_non_lesional <- clusterVolumes[,grep("H", colnames(cluster_and_samples)[c(1:(ncol(cluster_and_samples)-3))])]
  } else {
    cluster_and_samples$Volume_non_lesional <- rowSums(clusterVolumes[,grep("H", colnames(cluster_and_samples)[c(1:(ncol(cluster_and_samples)-3))])])
  }
  
  ## Loop over the genes and calculate the volume of the clusters that carry
  ## mutations in those genes as a fraction of the total volume for the patient.
  thisPat <- cell_frac_input[cell_frac_input$patient==patient & !cell_frac_input$isSubclone,]
  c <- 1
  for(gene in geneList$V1) {
    thisGene <- thisPat[thisPat$gene==gene,]
    if(dim(thisGene)[1] < 1) {
      # No mutations found in this gene in this patient. 
      cell_frac_results_total[r,c] <- 0
      cell_frac_results_lesional[r,c] <- 0
      cell_frac_results_non_lesional[r,c] <- 0
      
    } else {
      
      cell_frac_results_total[r,c] <- sum(cluster_and_samples$Volume[cluster_and_samples$cluster_id %in% thisGene$ClusterID])/sum(volumes$CutVolume)
      cell_frac_results_lesional[r,c] <- sum(cluster_and_samples$Volume_lesional[cluster_and_samples$cluster_id %in% thisGene$ClusterID])/sum(volumes$CutVolume[grep("H", volumes$SampleID, invert=T)])
      cell_frac_results_non_lesional[r,c] <- sum(cluster_and_samples$Volume_non_lesional[cluster_and_samples$cluster_id %in% thisGene$ClusterID])/sum(volumes$CutVolume[grep("H", volumes$SampleID)])
      
    }
    c <- c+1
  }
  r <- r+1
}

cell_frac_results_total <- as.data.frame(cell_frac_results_total)
cell_frac_results_lesional <- as.data.frame(cell_frac_results_lesional)
cell_frac_results_non_lesional <- as.data.frame(cell_frac_results_non_lesional)
colnames(cell_frac_results_total) = colnames(cell_frac_results_lesional) = colnames(cell_frac_results_non_lesional) <- geneList$V1
rownames(cell_frac_results_total) =rownames(cell_frac_results_lesional) = rownames(cell_frac_results_non_lesional) <- patient_list

## When no mutations are found for a patient, I give the value 0. This only gets set as NA when
## There are mutations but no non-lesional samples. 
cell_frac_results_non_lesional[!complete.cases(cell_frac_results_non_lesional),] <- NA

write.table(cell_frac_results_total, file="/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/mutant_cell_fractions.txt", sep="\t", row.names = T, col.names = T, quote = F)
write.table(cell_frac_results_lesional, file="/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/mutant_cell_fractions_lesional.txt", sep="\t", row.names = T, col.names = T, quote = F)
write.table(cell_frac_results_non_lesional, file="/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/mutant_cell_fractions_non_lesional.txt", sep="\t", row.names = T, col.names = T, quote = F)


library(ggplot2)
library(reshape2)
cell_frac_results_m <- melt(cell_frac_results_total)

ggplot(cell_frac_results_m, aes(x=value, fill=variable)) + geom_histogram(bins=50)+
  facet_wrap(.~variable, nrow=2) 





