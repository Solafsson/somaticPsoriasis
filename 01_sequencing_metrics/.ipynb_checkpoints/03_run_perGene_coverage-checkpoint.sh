#!/bin/bash

## The purpose of this script is to take in a bed file containing the coordinates of the coding regions
## of the genes found to be under positive selection and calculate the coverage of each sample over those 
## coordinates. This is so I can decide which samples have enough reads of the genes for it to make sense
## to include those samples in the estimates of the fraction of mutated cells. 


output_dir=/lustre/scratch119/humgen/projects/psoriasis/coverage_calculation/
script_dir=/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/01_sequencing_metrics/
clusterOutputDir=${output_dir}logs/
mkdir -p ${clusterOutputDir}

# Where the bam files and the CaveMan/Pindel outputs get stored
pso_dataDir=/nfs/cancer_ref01/nst_links/live/


metaData=$1
awk '($8=="WGS" || $8=="WES") && ($16=="NA" || $16=="Technical_duplicate") {print $1, $17}' < ${metaData} > ${output_dir}sample_list.tmp


# Use the coordinates of the genes according to the UCSC genome browser
# I want to intersect those coordinates with the bed file for the bait set so I'm only calculating
# the coverage over the regions that were captured in the bait-set. 
gene_list=${script_dir}coding_regions_genes_selection.bed

while read chr startPos endPos gene; do
    echo "${chr} ${startPos} ${endPos} ${gene}" | sed 's/ /\t/g' > ${output_dir}tmp
    /software/team152/bedtools intersect -a /lustre/scratch119/humgen/projects/psoriasis/resources/Human_all_exon_V5.bed -b ${output_dir}tmp > ${output_dir}${gene}_bait_regions.bed
    
    gene_coordinates=${output_dir}${gene}_bait_regions.bed
    
    while read sample projectNr; do
        bsub -o ${clusterOutputDir}/${sample}_samtools.out -e ${clusterOutputDir}/${sample}_samtools.err -R"select[mem>200] rusage[mem=200]" -M200 \
        "bash ${script_dir}perGene_coverage_worker.sh ${sample} ${projectNr} ${output_dir} ${gene_coordinates} ${gene}" 
    done < ${output_dir}sample_list.tmp

done < ${gene_list}

## join -1 1 -2 1 <( sort -k1,1 samtools_CHEK2_coverage.txt) <(sort -k1,1 samtools_EEF1A1_coverage.txt) | join -1 1 -2 1 - <(sort -k1,1 samtools_FAT1_coverage.txt) | join -1 1 -2 1 - <(sort -k1,1 samtools_GXYLT1_coverage.txt) | join -1 1 -2 1 - <(sort -k1,1 samtools_MRC1_coverage.txt) | join -1 1 -2 1 - <(sort -k1,1 samtools_NOTCH1_coverage.txt) | join -1 1 -2 1 - <(sort -k1,1 samtools_NOTCH2_coverage.txt) | join -1 1 -2 1 - <(sort -k1,1 samtools_PPM1D_coverage.txt) | join -1 1 -2 1 - <(sort -k1,1 samtools_TP53_coverage.txt) | join -1 1 -2 1 - <(sort -k1,1 samtools_ZFP36L2_coverage.txt) | awk 'BEGIN {print "sampleID", "CHEK2_cov", "EEF1A1_cov", "FAT1_cov", "GXYLT1_cov", "MRC1_cov", "NOTCH1_cov", "NOTCH2_cov", "PPM1D_cov", "TP53_cov", "ZFP36L2_cov"} {print $1, $7, $13, $19, $25, $31, $37, $43, $49, $55, $61}' > gene_coverage.txt



exit $?






