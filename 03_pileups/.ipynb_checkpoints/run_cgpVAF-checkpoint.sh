#!/bin/bash

## See this link for instructions on how to run cgpVAF:
## https://confluence.sanger.ac.uk/pages/viewpage.action?pageId=22710418
##        Run /software/CGP/projects/cgpVAFcommand/perl/bin/createVafCmd.pl -help for documentation

## Usage: cgpVAF is a hellish software that fails when it god damn wants to and doesn't need to give any reason.
##        This script runs reasonably reliably.
##        run_cgpVAF.sh <patient> <project_nr> <mutationType>
##        Example:
##        for i in {7..19}; do  ~/phd/somatic_ibd_p1/pileups/run_cgpVAF.sh patient${i} 1874 snp; done

## Input: patient is a patient ID found in the sample_meta file. Project nr is the CanApps ID for the project.
##        This is 1874 for the IBD project and 1494 or 1728 for the normal data.
##        mutationType is either snp or indel. If mutationType=snp, then I expect this file to exist:
##       ${hairpin_filter_dir}${patient}/${crypt}_complete_final_retained_3.vcf
##       (Change that so I can take the file outputted by the hairpin filter)

## /nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/03_pileups/run_cgpVAF.sh patient01 2545 indel

## If cgpVAF refuses to acknowledge the -be flag then you can just rename the bam and .bai files in the bam file directory
## to be just <sample>.bam
## rename 's/.sample.dupmarked//' *.bam
## rename 's/.sample.dupmarked//' *.bam.bai

## samtools view -C <input.bam> -T <reference.fasta> -o <output.cram>


output_dir=/lustre/scratch119/humgen/projects/psoriasis/pileups/
script_dir=/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/03_pileups/

hairpin_filter_dir=/lustre/scratch119/humgen/projects/psoriasis/hairpin_filters/
#reference=/lustre/scratch119/casm/team78pipelines/canpipe/live/ref/human/GRCh37d5/genome.fa
reference=/lustre/scratch112/sanger/cgppipe/canpipe/test/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa
clusterOutputDir=${output_dir}logs/

patient=$1
mutationType=$2
sample_meta=$3

pso_dataDir=/nfs/cancer_ref01/nst_links/live/

mkdir -p ${clusterOutputDir}
mkdir -p ${output_dir}
mkdir -p ${output_dir}${patient}
mkdir -p ${output_dir}bam_file_inputs
mkdir -p ${output_dir}bed_file_inputs
mkdir -p ${output_dir}bam_file_inputs/${patient}
mkdir -p ${output_dir}bed_file_inputs/${patient}

rm -f ${output_dir}${patient}/${patient}_err
rm -f ${output_dir}${patient}/${patient}_log

awk -v PAT="$patient" '$6==PAT && ($8=="WGS" || $8=="WES") && $12!="Matched_Normal" && ($16=="NA" || $16=="Technical_Duplicate") {print $1, $17}' < ${sample_meta} > ${output_dir}${patient}/${patient}_sampleList.txt
awk -v PAT="$patient" '$6==PAT && ($8=="WGS" || $8=="WES") && $12!="Matched_Normal" && ($16=="NA" || $16=="Technical_Duplicate") && $15!="NA" {print $15, $17}' < ${sample_meta} | sort -u > ${output_dir}${patient}/${patient}_normalList.txt

while read norm project_nr; do
    ## Make links to normal bam files
    ln -s ${pso_dataDir}${project_nr}/${norm}/${norm}.sample.dupmarked.bam ${output_dir}bam_file_inputs/${patient}/
    ln -s ${pso_dataDir}${project_nr}/${norm}/${norm}.sample.dupmarked.bam.bai ${output_dir}bam_file_inputs/${patient}/
done < ${output_dir}${patient}/${patient}_normalList.txt

rm -f ${output_dir}bed_file_inputs/${patient}_indels.bed

while read crypt project_nr; do
        ln -s ${pso_dataDir}${project_nr}/${crypt}/${crypt}.sample.dupmarked.bam ${output_dir}bam_file_inputs/${patient}/
        ln -s ${pso_dataDir}${project_nr}/${crypt}/${crypt}.sample.dupmarked.bam.bai ${output_dir}bam_file_inputs/${patient}/


        ## Make bed file although not actually .bed - chr pos ref alt
        if [ "$mutationType" == "snp" ]; then
            grep -v -e '^\#' ${hairpin_filter_dir}${patient}/${crypt}_complete_final_retained_3.vcf | awk '{OFS="\t"; print $1, $2, $4, $5}' \
            > ${output_dir}bed_file_inputs/${patient}/${crypt}_filteredMutations.bed

        elif [ "$mutationType" == "indel" ]; then
        ## Filter the homopolimers
        #perl ~/phd/somatic_ibd_p1/qc/filter_homopolymers/indel_homopolymer_filter.pl <( zcat ${pso_dataDir}${project_nr}/${crypt}/${crypt}.pindel.annot.vcf.gz ) | \
         zcat ${pso_dataDir}${project_nr}/${crypt}/${crypt}.pindel.annot.vcf.gz | \
         grep -v -e '^\#' | awk '{print $1":"$2":"$4":"$5":"$11}' \
            | awk 'BEGIN {FS=":"} $NF!="." {OFS="\t"; print $1, $2, $3, $4}' >> ${output_dir}bed_file_inputs/${patient}_indels.bed
        fi
done < ${output_dir}${patient}/${patient}_sampleList.txt


# Combine and sort bed files
if [ "$mutationType" == "snp" ]; then
    cat ${output_dir}bed_file_inputs/${patient}/*_filteredMutations.bed | sort -u > ${output_dir}bed_file_inputs/${patient}_filteredMutations.bed
    rm -r ${output_dir}bed_file_inputs/${patient}/
elif [ "$mutationType" == "indel" ]; then
    awk '$1=="Y" || $1=="X" {print}' < ${output_dir}bed_file_inputs/${patient}_indels.bed > ${output_dir}bed_file_inputs/${patient}_sexChr.tmp
    awk '$1!="Y" && $1!="X" {print}' < ${output_dir}bed_file_inputs/${patient}_indels.bed > ${output_dir}bed_file_inputs/${patient}_autosomes.tmp
    sort -u ${output_dir}bed_file_inputs/${patient}_autosomes.tmp  | sort -n -k1,2 > ${output_dir}bed_file_inputs/${patient}_indels.bed
    sort -u ${output_dir}bed_file_inputs/${patient}_sexChr.tmp | sort -n -k1,2  >>  ${output_dir}bed_file_inputs/${patient}_indels.bed
    rm ${output_dir}bed_file_inputs/${patient}_autosomes.tmp
    rm ${output_dir}bed_file_inputs/${patient}_sexChr.tmp
fi


## If there isn't a normal then just run it with one of the samples as normal. cgpVaf doesn't actually require
## a matched normal.
nrLines=$( wc -l ${output_dir}${patient}/${patient}_normalList.txt | awk '{print $1}' )
if [ "$nrLines" -gt 0 ]; then
    normal=$( head -1 ${output_dir}${patient}/${patient}_normalList.txt | awk '{print $1}' )
else
    normal=$( head -1 ${output_dir}${patient}/${patient}_sampleList.txt | awk '{print $1}' )
    tail -n+2 ${output_dir}${patient}/${patient}_sampleList.txt > ${output_dir}${patient}/${patient}_tmp
    mv ${output_dir}${patient}/${patient}_tmp ${output_dir}${patient}/${patient}_sampleList.txt
fi

sampleList=$( awk '{ORS=","; print $1}' < ${output_dir}${patient}/${patient}_sampleList.txt | awk '{print substr($1, 1, length($1)-1)}' )


if [ "$mutationType" == "indel" ]; then
    mkdir -p ${output_dir}${patient}/indel/
    echo "bsub -q basement -o ${output_dir}${patient}/${patient}_log -e ${output_dir}${patient}/${patient}_err \
    'bash ${script_dir}cgpVaf_worker.sh ${output_dir}bam_file_inputs/${patient}/ ${reference} ${output_dir}${patient}/indel/ ${normal} ${sampleList} ${output_dir}bed_file_inputs/${patient}_indels.bed ${mutationType}'"

    bsub -q basement -o ${output_dir}${patient}/${patient}_log -e ${output_dir}${patient}/${patient}_err -R"select[mem>2000] rusage[mem=2000]" -M2000 \
    "bash ${script_dir}cgpVaf_worker.sh ${output_dir}bam_file_inputs/${patient}/ ${reference} ${output_dir}${patient}/indel/ ${normal} ${sampleList} ${output_dir}bed_file_inputs/${patient}_indels.bed ${mutationType}"
elif [ "$mutationType" == "snp" ]; then
    mkdir -p ${output_dir}${patient}/snp/
        echo "bsub -o ${output_dir}${patient}/${patient}_log -e ${output_dir}${patient}/${patient}_err \
    'bash ${script_dir}cgpVaf_worker.sh ${output_dir}bam_file_inputs/${patient}/ ${reference} ${output_dir}${patient}/snp/ ${normal} ${sampleList} ${output_dir}bed_file_inputs/${patient}_filteredMutations.bed ${mutationType}'"

    bsub -q long -o ${output_dir}${patient}/${patient}_log -e ${output_dir}${patient}/${patient}_err -R"select[mem>2000] rusage[mem=2000]" -M2000 \
    "bash ${script_dir}cgpVaf_worker.sh ${output_dir}bam_file_inputs/${patient}/ ${reference} ${output_dir}${patient}/snp/ ${normal} ${sampleList} ${output_dir}bed_file_inputs/${patient}_filteredMutations.bed ${mutationType}"
fi


exit $?