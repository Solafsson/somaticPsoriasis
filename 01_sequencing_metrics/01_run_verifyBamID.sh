#!/bin/bash

output_dir=/lustre/scratch119/humgen/projects/psoriasis/contamination_check/
#metaData=/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/sample_info/metaData.txt
#metaData=${output_dir}contamination_sample_list_01102019.txt
ThousandGenomes_variants=/lustre/scratch119/humgen/projects/psoriasis/resources/
#verifyBamID_exec=/lustre/scratch116/casm/cgp/users/fa8/bin/verifyBamID
verifyBamID_exec=/lustre/scratch119/casm/team268im/fa8/VerifyBamId-farm5/VerifyBamID
bamDir=/nfs/cancer_ref01/nst_links/live/
script_dir=/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/01_sequencing_metrics/
REFERENCE=/lustre/scratch117/core/sciops_repository/references/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/all/fasta/Homo_sapiens.GRCh38_full_analysis_set_plus_decoy_hla.fa
sample_list=$1
#sample_list should contain two columns. The sampleID and the CanApps project number of the sample.
#awk 'NR>1 {print $1, $17}' < ${metaData} > ${output_dir}sample_list.txt

## Now using version 2 of VerifyBamID:
## https://github.com/Griffan/VerifyBamID
## Reference this paper:
## https://genome.cshlp.org/content/30/2/185


nrJobs=$( wc -l ${sample_list} | awk '{print $1}' )

mkdir -p ${output_dir}logs

bsub -J"check_contamination[1-${nrJobs}]" -q long -M3000 -R'span[hosts=1] select[mem>3000] rusage[mem=3000]' \
-e ${output_dir}logs/check_contamination_errors.%J.%I -o ${output_dir}logs/check_contamination_output.%J.%I \
bash ${script_dir}verifyBamID_worker.sh ${sample_list} ${output_dir} ${ThousandGenomes_variants} \
 ${verifyBamID_exec} ${bamDir} ${REFERENCE}


exit $?