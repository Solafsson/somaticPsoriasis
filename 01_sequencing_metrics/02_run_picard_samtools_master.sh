#!/bin/bash

output_dir=/lustre/scratch119/humgen/projects/psoriasis/coverage_calculation/hsmetrics/
script_dir=/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/01_sequencing_metrics/
clusterOutputDir=${output_dir}logs/
mkdir -p ${clusterOutputDir}

# Where the bam files and the CaveMan/Pindel outputs get stored
pso_dataDir=/nfs/cancer_ref01/nst_links/live/

metaData=$1


awk '($8=="WGS" || $8=="WES") && ($16=="NA" || $1=="Technical_Duplicate") {print $1, $17}' < ${metaData} > ${output_dir}sample_list.tmp


while read sample projectNr; do

   # bsub -o ${clusterOutputDir}/${sample}_samtools.out -e ${clusterOutputDir}/${sample}_samtools.err -R"select[mem>200] rusage[mem=200]" -M200 \
  #  "bash ${script_dir}samtools_depth_worker.sh ${sample} ${projectNr} ${output_dir}" 

    bsub -o ${clusterOutputDir}/${sample}_picard.out -e ${clusterOutputDir}/${sample}_picard.err -R"select[mem>6000] rusage[mem=6000]" -M6000 \
    "bash ${script_dir}picard_metrics_worker.sh ${sample} ${projectNr} ${output_dir}" 

done < ${output_dir}sample_list.tmp


exit $?