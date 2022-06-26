#!/bin/bash

## The purpose of this script is to handle the submission of the pileup work to the farm. 
## /nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/03_pileups/run_pileup_bam2R.sh  /nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/sample_info/sample_meta.txt

## Input: metaData is a sample-metaData file and should have the following format:
#sampleID	Date_of_dissection	SANGER_PLATE_ID	SANGER_SAMPLE_ID	biopsy_ID	patient_ID	Library_concentration	Seq_status	Volume	Area	Location	type	lcm_date	Coverage	pipeline_normal	excl_criteria	projectNr
#P01H_LL_6	05.02.2021	DN774429I	6308STDY9901922	P01H	patient1	70.409	WES	661	9100	NA	Healthy	05.02.2021	NA	P04H_T_6	NA	2545
#P01H_LL_7	05.02.2021	DN774429I	6308STDY9901923	P01H	patient1	51.14	Wait	272	4700	NA	Healthy	05.02.2021	NA	P04H_T_6	NA	2545
#P01H_LL_8	05.02.2021	DN774429I	6308STDY9901924	P01H	patient1	24.431	Wait	579	6600	NA	Healthy	05.02.2021	NA	P04H_T_6	NA	2545


#           patient_subset is a comma-separated list of patient_IDs, matching values found in the sample-metaData file. Analysis will be restricted
#           to those patients. If nothing is specified, analysis is run on all unique patient values in the metaData file.
#           output_dir is a directory on lustre. If nothing is specified results are written to /lustre/scratch119/humgen/projects/psoriasis/pileups/.
#           clusterErrorOutputDir is a directory to which the log files of the jobs will be written.


output_dir=/lustre/scratch119/humgen/projects/psoriasis/pileups/
script_dir=/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/03_pileups/
clusterOutputDir=${output_dir}logs/
mkdir -p ${clusterOutputDir}

# Where the bam files and the CaveMan/Pindel outputs get stored
pso_dataDir=/nfs/cancer_ref01/nst_links/live/

metaData=$1

# An optional second argument. patient_subset is a comma-separated list of patient_IDs, matching values found in the sample-metaData file. Analysis will be restricted  to those patients. If nothing is specified, analysis is run on all unique patient values in the metaData file.
patient_subset=$2


# Create a list of subjects
# Restrict to patient_subset if applicable. 
awk '($8=="WGS" || $8=="WES") && ($16=="NA" || $1=="Technical_Duplicate") {print $6}' < ${metaData} | sort -u > ${output_dir}patientList_pileup.tmp
if [ ! -z "$patient_subset" ]; then
    echo ${patient_subset} | awk -F, -v OFS="\n" '{$1=$1; print}' > ${output_dir}patient_subset.tmp
    grep -w -f ${output_dir}patient_subset.tmp ${output_dir}patientList_pileup.tmp > ${output_dir}tmp
    mv ${output_dir}tmp ${output_dir}patientList_pileup.tmp
fi


## Create comma-separated vectors of project numbers and the corresponding samples
## Then submit one job to the farm for each patient:
while read patient; do
    mkdir -p ${output_dir}${patient}/
    projectVector=$( awk -v PAT="$patient" 'BEGIN {projectV=""} $6==PAT && ($8=="WGS" || $8=="WES") && ($16=="NA" || $16=="Technical_Duplicate") {projectV=projectV","$17} END {print substr(projectV, 2, length(projectV))}' < ${metaData} )

    sampleVector=$( awk -v PAT="$patient" 'BEGIN {sampleV=""} $6==PAT && ($8=="WGS" || $8=="WES") && ($16=="NA" || $16=="Technical_Duplicate") {if(sampleV=="") {sampleV=$1} else {sampleV=sampleV","$1}} END {print sampleV}' < ${metaData} )

echo "Submitting to farm:"
echo "bsub -o ${clusterOutputDir}${patient}.out -e ${clusterOutputDir}${patient}.err -q long -R'select[mem>3000] rusage[mem=3000]" -M3000 "/software/R-4.1.0/bin/Rscript ${script_dir}pileup_worker_bam2R.R' ${patient} ${sampleVector} ${projectVector} ${output_dir} ${pso_dataDir}"    
    bsub -o ${clusterOutputDir}${patient}.out -e ${clusterOutputDir}${patient}.err -q long -R"select[mem>3000] rusage[mem=3000]" -M3000 "/software/R-4.1.0/bin/Rscript ${script_dir}pileup_worker_bam2R_indels.R" ${patient} ${sampleVector} ${projectVector} ${output_dir} ${pso_dataDir}

done < ${output_dir}patientList_pileup.tmp


exit $?