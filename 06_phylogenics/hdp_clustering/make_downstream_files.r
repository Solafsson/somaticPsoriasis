
hdp_dir="/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/"
patient_list <- read.table(paste(hdp_dir,"patientList_hdp.tmp", sep=""))

# failed: 6,7,110,55,79
for(patient in patient_list$V1[93:112]) {
  load(paste(hdp_dir, patient, "/", "ndp_", patient, "_2021_12_09/Rsession_ls2.dat", sep=""))
  
  mut_assignment <- data.frame(Chr=ls_tbl$chrom, Pos=ls_tbl$pos, Ref=lymph.depth$ref,
                               Alt=lymph.depth$alt, Cluster=ls_tbl$clust_assign,Patient=patient,
                               ClusterID=paste(patient, ls_tbl$clust_assign, sep="_"))
  
  write.table(mut_assignment, paste(hdp_dir, patient, "/", "ndp_", patient, "_2021_12_09/sbs_assignment.txt", sep=""), sep="\t", row.names = F, quote=F)

  
}




