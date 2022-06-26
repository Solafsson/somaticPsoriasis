#!/bin/bash

## The purpose of this script is to handle the job submission for 011_indel_posterior_sampling_worker.r to the farm.


mut_type="indels"
nr_chains=20

script_dir=/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/07_signature_extraction/
output_dir=/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/hdp/${mut_type}/
clusterOutputDir=${output_dir}logs/
mkdir -p ${clusterOutputDir}


for i in $(seq 1 $nr_chains); do

    bsub -o ${clusterOutputDir}hdp_${i}.out -e ${clusterOutputDir}hdp_${i}.err -q long -n4 -R 'select[mem>=3000] rusage[mem=3000]' -M3000 \
  "/software/R-4.1.0/bin/Rscript ${script_dir}011_indel_posterior_sampling_worker.r ${mut_type} ${output_dir} ${i}" 
done



exit $?