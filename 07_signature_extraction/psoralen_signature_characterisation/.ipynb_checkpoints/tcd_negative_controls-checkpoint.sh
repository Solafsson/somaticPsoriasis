#!/bin/bash

## The purpose of this script is to repeat the analyses documented in puva_tcd.sh but for
## other mutation types. This is to compare the transcription-coupled damage we see for PUVA with 
## That of other mutation types, including those characterizing UV-exposure and C>A mutations
## (which are a feature of neither.)

## I will look at C>T (or G>A) mutations and compare these with any C>A (or G>T) mutation.

# Compile a list of genes highly expressed in the skin. Include strand information.
# Create 1kb bins 10kb up and downstream of the TSSs
# Count the number of C>T and C>A mutations on the transcribed strand
# and the number of G>A and G>T mutations on the untranscribed strand
# if the gene is on the (+) strand, the transcribed strand is the complement of the 
# reference strand but if the gene is on the (-) strand, the transcribed strand is the reference. 

# That means:
# To find the number of mutations overlapping the transcribed strand:
# If the gene is on the (+) strand, count the number of G>A and G>T mutations overlapping it.
# If the gene is on the (-) strand, count the C>T and C>A mutations. 

# To find the number of mutations overlapping the non-transcribed strand, do the reverse: 
# If the gene is on the (+) strand, count the C>T and C>A mutations 

 

## I downloaded the GTEx expression data: GTEx Analysis V8 (dbGaP Accession phs000424.v8.p2)
## from https://www.gtexportal.org/home/datasets in November 2021

#bsub -o tcd.out -e tcd.err -q normal -R"select[mem>5000] rusage[mem=5000]" -M5000 "bash  ${pso_nfs}07_signature_extraction/strand_asymmetries/tcd_negative_controls.sh"

bedtools_exec=/software/team152/bedtools2/bin/bedtools
puva_dir=/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/
reference=/lustre/scratch117/core/sciops_repository/references/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/all/fasta/Homo_sapiens.GRCh38_full_analysis_set_plus_decoy_hla.fa


## Count the number of AT/TA dinucleotide sites in each bin:
# May do this later to standardize the mutation rate by sequence content but not at the moment. 
#${bedtools_exec} nuc -fi ${reference} -bed ${puva_dir}top10pc_expr_genes_skin.bed \
#   -pattern AT -C > ${puva_dir}top10pc_expr_genes_skin_AT_counts.txt

#${bedtools_exec} nuc -fi ${reference} -bed ${puva_dir}top10pc_expr_genes_skin.bed \
#    -pattern TA -C > ${puva_dir}top10pc_expr_genes_skin_TA_counts.txt



rm -f ${puva_dir}*_Guv_transcribed.bed
rm -f ${puva_dir}*_Cuv_transcribed.bed
rm -f ${puva_dir}*_Gother_transcribed.bed
rm -f ${puva_dir}*_Cother_transcribed.bed
while read sample; do

    awk '$4=="C" && $5=="T" {OFS="\t"; print $1, $2, $3, $4, $5}' < ${puva_dir}${sample}_vcf_no_header.txt > ${puva_dir}${sample}_UV_mutations_refC.bed
    
    awk '$4=="G" && $5=="A" {OFS="\t"; print $1, $2, $3, $4, $5}' < ${puva_dir}${sample}_vcf_no_header.txt > ${puva_dir}${sample}_UV_mutations_refG.bed
    
    awk '$4=="C" && ($5=="A" || $5=="G") {OFS="\t"; print $1, $2, $3, $4, $5}' < ${puva_dir}${sample}_vcf_no_header.txt > ${puva_dir}${sample}_other_mutations_refC.bed
    
    awk '$4=="G" && ($5=="T" || $5=="C") {OFS="\t"; print $1, $2, $3, $4, $5}' < ${puva_dir}${sample}_vcf_no_header.txt > ${puva_dir}${sample}_other_mutations_refG.bed


    ## Find the mutations occurring on the transcribed strand. 
    ## First G>A and G>T mutations in genes falling on the positive strand
    ${bedtools_exec} intersect -a ${puva_dir}top10pc_expr_genes_skin_posStrand.bed -b ${puva_dir}${sample}_UV_mutations_refG.bed >> ${puva_dir}${sample}_Guv_transcribed.bed
    ${bedtools_exec} intersect -a ${puva_dir}top10pc_expr_genes_skin_posStrand.bed -b ${puva_dir}${sample}_other_mutations_refG.bed >> ${puva_dir}${sample}_Gother_transcribed.bed
    
    # And next C>T and C>A mutations in genes falling on the negative strand
    ${bedtools_exec} intersect -a ${puva_dir}top10pc_expr_genes_skin_negStrand.bed -b ${puva_dir}${sample}_UV_mutations_refC.bed >> ${puva_dir}${sample}_Cuv_transcribed.bed
    ${bedtools_exec} intersect -a ${puva_dir}top10pc_expr_genes_skin_negStrand.bed -b ${puva_dir}${sample}_other_mutations_refC.bed >> ${puva_dir}${sample}_Cother_transcribed.bed

    ## Next find the mutations occurring on the non-transcribed strand.
    ## C>T and C>A mutations on the positive strand
    ${bedtools_exec} intersect -a ${puva_dir}top10pc_expr_genes_skin_posStrand.bed -b ${puva_dir}${sample}_UV_mutations_refC.bed >> ${puva_dir}${sample}_Cuv_transcribed.bed
    ${bedtools_exec} intersect -a ${puva_dir}top10pc_expr_genes_skin_posStrand.bed -b ${puva_dir}${sample}_other_mutations_refC.bed >> ${puva_dir}${sample}_Cother_transcribed.bed

    ${bedtools_exec} intersect -a ${puva_dir}top10pc_expr_genes_skin_negStrand.bed -b ${puva_dir}${sample}_UV_mutations_refG.bed >> ${puva_dir}${sample}_Guv_transcribed.bed
    ${bedtools_exec} intersect -a ${puva_dir}top10pc_expr_genes_skin_negStrand.bed -b ${puva_dir}${sample}_other_mutations_refG.bed >> ${puva_dir}${sample}_Gother_transcribed.bed

done < ${puva_dir}sample_list.txt

cat ${puva_dir}*_Cuv_transcribed.bed > ${puva_dir}combined_Cuv_transcribed.bed
cat ${puva_dir}*_Guv_transcribed.bed > ${puva_dir}combined_Guv_transcribed.bed
cat ${puva_dir}*_Cother_transcribed.bed > ${puva_dir}combined_Cother_transcribed.bed
cat ${puva_dir}*_Gother_transcribed.bed > ${puva_dir}combined_Gother_transcribed.bed

exit $?