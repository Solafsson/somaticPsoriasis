
## The purpose of this script is to plot fraction of shared mutations against
## distance and to compare the extent of long-distance clone expansion between
## patients who show or don't show the PUVA mutational signature.




### Read in all pairwise distances
spatial_mat_dir="/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/sample_info/spatial_relationships/"
binomial_dir="/lustre/scratch119/humgen/projects/psoriasis/binomial_filters/"
files <- list.files(spatial_mat_dir)
patients <- as.character()
frac_shared <- as.numeric()
distances <- as.numeric()
pairs <- as.character()
k <- 1
for(f in files) {

  patient <- gsub("P", "patient", unlist(strsplit(f, split="H|L"))[c(T,F)])
  tmp <- read.table(paste(spatial_mat_dir, f, sep=""), h=T)
  
  genotype_bin <- read.table(paste(binomial_dir, patient, "/", patient, "_genotype_allMuts.txt", sep=""), h=T)
  for(i in 1:ncol(tmp)) {
    for(j in 1:nrow(tmp)) {
      if(!is.na(tmp[j,i]) & tmp[j,i]>0) {
        patients[k] <- patient
        frac_shared[k] <- sum(genotype_bin[,colnames(tmp)[i]]>0 & genotype_bin[,rownames(tmp)[j]]>0)/sum(genotype_bin[,colnames(tmp)[i]]>0 | genotype_bin[,rownames(tmp)[j]]>0)
        distances[k] <- tmp[j,i]
        pairs[k] <- paste(colnames(tmp)[i], rownames(tmp)[j], sep=",")
        
        k <- k+1
      }
    }
  }
}

shared_by_distance <- data.frame(patients, pairs, frac_shared, distances)
shared_by_distance$PUVA_exposed <- F
shared_by_distance$PUVA_exposed[shared_by_distance$patients %in% puva_exp_patients] <- T

write.table(shared_by_distance, file="/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/shared_by_distance.txt", sep="\t", row.names = F, quote=F)


