#!/bin/bash

## The purpose of this script is to test if there is an effect of replication timing on 
## PUVA-associated mutagenesis. 
## I have got four bed files, corresponding to bins of replication timing. I want to estimate
## the per-base mutation rate in each 

## This is going to have the following steps:
## 1. Create a .bed file from my mutation calls. Only include T>A, T>C and T>G mutations. 
## 2. Count the number of TA/AT dinucleotides in each genomic segment in each .bed file
## 3. Calculate the number of mutations per TA/AT dinucleotides. 

## https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html
## https://bedtools.readthedocs.io/en/latest/content/tools/nuc.html
## https://www.biostars.org/p/184842/

#bsub -o RT.out -e RT.err -q normal -R"select[mem>5000] rusage[mem=5000]" -M5000 "bash  ${pso_nfs}07_signature_extraction/strand_asymmetries/puva_replication_timing.sh"

bedtools_exec=/software/team152/bedtools2/bin/bedtools
puva_dir=/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/
reference=/lustre/scratch117/core/sciops_repository/references/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/all/fasta/Homo_sapiens.GRCh38_full_analysis_set_plus_decoy_hla.fa
bedFile_dir=/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/07_signature_extraction/strand_asymmetries/

    ## 2. Count the number of TA/AT dinucleotides in each genomic segment in each .bed file
for i in {1..4}; do

    if [ ! -f ${puva_dir}RT_${i}_hg38_AT_and_TA_counts.txt ]; then
        ${bedtools_exec} nuc -fi ${reference} -bed ${bedFile_dir}RT_${i}_hg38.bed \
        -pattern AT -C > ${puva_dir}RT_${i}_hg38_AT_counts.txt

        ${bedtools_exec} nuc -fi ${reference} -bed ${bedFile_dir}RT_${i}_hg38.bed \
        -pattern TA -C > ${puva_dir}RT_${i}_hg38_TA_counts.txt

        paste ${puva_dir}RT_${i}_hg38_AT_counts.txt ${puva_dir}RT_${i}_hg38_TA_counts.txt \
        | awk -F "\t" 'BEGIN {OFS="\t"; print "Chr", "Start", "End", "AT_count", "TA_count"} NR>1 {OFS="\t"; print $1, $2, $3,$13,$26}' \
        > ${puva_dir}RT_${i}_hg38_AT_and_TA_counts.txt
        
        rm ${puva_dir}RT_${i}_hg38_AT_counts.txt
        rm ${puva_dir}RT_${i}_hg38_TA_counts.txt
        
    fi
done

while read sample; do

    ## 1. Create a .bed file from my mutation calls. 
    grep -v '\#' ${puva_dir}${sample}.vcf | awk -F "\t" '$4=="A" || $4=="T" {OFS="\t"; print $1, $2-1, $3-1, $4, $5}' > ${puva_dir}${sample}_PUVA_mutations.bed

    rm -f ${puva_dir}${sample}_RT_results.txt
    for i in {1..4}; do
        nrMuts=$( ${bedtools_exec} intersect -a ${bedFile_dir}RT_${i}_hg38.bed -b ${puva_dir}${sample}_PUVA_mutations.bed | wc -l )
        mutRate=$( awk -v muts="$nrMuts" 'BEGIN {sum=0} NR>1 {sum=sum+$4+$5} END {print muts/sum}' < ${puva_dir}RT_${i}_hg38_AT_and_TA_counts.txt )
        echo "RT_${i} ${mutRate}" >> ${puva_dir}${sample}_RT_results.txt

    done
done < ${puva_dir}sample_list.txt

exit $?