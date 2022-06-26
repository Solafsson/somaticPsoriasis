#!/bin/bash

## The purpose of this script is to take in a sample-metadata file and run the beta-binomial filters for all samples
## in the file.
## The script will also plot heatmaps of the VAFs, the VAF distribution of each sample, write a binary genotype matrix for the patient and write out the NV, NR and VAFs of all passing mutations for each patient. 


## run_binomial_filters sample_metaData patient_metaData patient_subset
## Example:
## /nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/04_binomial_filters/run_binomial_filters.sh /nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/sample_info/sample_meta.txt /nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/sample_info/patient_meta.txt patient01


sample_metaData=$1
patient_metaData=$2
patient_subset=$3

## some hard-coded variables:
pileup_dir=/lustre/scratch119/humgen/projects/psoriasis/pileups/
script_dir=/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/04_binomial_filters/
output_dir=/lustre/scratch119/humgen/projects/psoriasis/binomial_filters/
clusterErrorOutput=${output_dir}logs/

mkdir -p ${clusterErrorOutput}


# Create a list of subjects
awk '( $8=="WGS" || $8=="WES" ) && ($16=="NA" || $16=="Technical_Duplicate") {print $6}' < ${sample_metaData} | sort -u > ${output_dir}patientList_binomialFiltering.tmp
if [ ! -z "$patient_subset" ]; then
    echo ${patient_subset} | awk -F, -v OFS="\n" '{$1=$1; print}' > ${output_dir}patient_subset.tmp
    grep -w -f ${output_dir}patient_subset.tmp ${output_dir}patientList_binomialFiltering.tmp > ${output_dir}tmp
    mv ${output_dir}tmp ${output_dir}patientList_binomialFiltering.tmp
fi


while read patient; do

    mkdir -p ${output_dir}${patient}

    ## What is the sex of the patient?
    sex=$( grep -w ${patient} ${patient_metaData} | awk '{print $5}' )
    if [ "$sex" == "male" ]; then
        sex=male
    else
        sex=female
    fi


    ## Submit R-script to the cluster

    bsub -q normal -o ${clusterErrorOutput}/${patient}_binomial_log -e ${clusterErrorOutput}/${patient}_binomial_err -R"select[mem>2000] rusage[mem=2000]" -M2000 \
    "/software/R-4.1.0/bin/Rscript ${script_dir}binomial_filters_worker.R ${patient} ${sex} ${script_dir} ${output_dir}${patient}/ ${pileup_dir}"



done < ${output_dir}patientList_binomialFiltering.tmp



exit $?