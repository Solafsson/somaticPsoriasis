#!/bin/bash

## The purpose of this script is to look for evidence of transcription coupled damage
## in samples that show the PUVA signature. 
## See figure 7A from http://www.cell.com/abstract/S0092-8674(15)01714-6

# Compile a list of genes highly expressed in the skin. Include strand information.
# Create 1kb bins 10kb up and downstream of the TSSs
# Count the number of T>A, T>C and T>G mutations at TpA sites on the transcribed strand
# and the number of A>T, A>G and A>C mutations at ApT sites on the untranscribed strand
# if the gene is on the (+) strand, the transcribed strand is the complement of the 
# reference strand but if the gene is on the (-) strand, the transcribed strand is the reference. 

# That means:
# To find the number of T>N or A>N mutations overlapping the transcribed strand:
# If the gene is on the (+) strand, count the number of A>T, A>G and A>C mutations overlapping it.
# If the gene is on the (-) strand, count the T>A, T>C and T>G mutations mutations. 

# To find the number of mutations overlapping the non-transcribed strand, do the reverse: 
# If the gene is on the (+) strand, count the T>A, T>C and T>G mutations mutations

## Do the same for other mutation types. Look at C>T (or G>A) mutations at
## TCC and CCC sites. These are the two top peaks of SBS7b. Compare with any C>A mutation. 

## I downloaded the GTEx expression data: GTEx Analysis V8 (dbGaP Accession phs000424.v8.p2)
## from https://www.gtexportal.org/home/datasets in November 2021

#bsub -o tcd.out -e tcd.err -q normal -R"select[mem>1000] rusage[mem=1000]" -M1000 "bash  ${pso_nfs}07_signature_extraction/strand_asymmetries/puva_tcd.sh"

bedtools_exec=/software/team152/bedtools2/bin/bedtools
puva_dir=/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/
reference=/lustre/scratch117/core/sciops_repository/references/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/all/fasta/Homo_sapiens.GRCh38_full_analysis_set_plus_decoy_hla.fa


## Count the number of AT/TA dinucleotide sites in each bin:
# May do this later to standardize the mutation rate by sequence content but not at the moment. 
#${bedtools_exec} nuc -fi ${reference} -bed ${puva_dir}top10pc_expr_genes_skin.bed \
#   -pattern AT -C > ${puva_dir}top10pc_expr_genes_skin_AT_counts.txt

#${bedtools_exec} nuc -fi ${reference} -bed ${puva_dir}top10pc_expr_genes_skin.bed \
#    -pattern TA -C > ${puva_dir}top10pc_expr_genes_skin_TA_counts.txt

for q in {1..5}; do

    ## This really only needed doing once:
    #if [ ! -f ${puva_dir}quintile${q}_expr_genes_skin.bed ]; then
        awk '$5=="+" {OFS="\t"; print $1, $2, $3, $4, $5, $6}' < ${puva_dir}quintile${q}_expr_genes_skin.bed > ${puva_dir}quintile${q}_expr_genes_skin_posStrand.bed

        awk '$5=="-" {OFS="\t"; print $1, $2, $3, $4, $5, $6}' < ${puva_dir}quintile${q}_expr_genes_skin.bed > ${puva_dir}quintile${q}_expr_genes_skin_negStrand.bed
    #fi 


    rm -f ${puva_dir}*_Tmut_transcribed.bed
    rm -f ${puva_dir}*_Amut_transcribed.bed
    while read sample; do

        awk '$4=="T" {OFS="\t"; print $1, $2, $3, $4, $5}' < ${puva_dir}${sample}_PUVA_mutations.bed > ${puva_dir}${sample}_PUVA_mutations_refT.bed

        awk '$4=="A" {OFS="\t"; print $1, $2, $3, $4, $5}' < ${puva_dir}${sample}_PUVA_mutations.bed > ${puva_dir}${sample}_PUVA_mutations_refA.bed

        ## Find the mutations occurring on the transcribed strand. 
        ## First A>N mutation in genes falling on the positive strand
        ${bedtools_exec} intersect -a ${puva_dir}quintile${q}_expr_genes_skin_posStrand.bed -b ${puva_dir}${sample}_PUVA_mutations_refA.bed >> ${puva_dir}${sample}_Tmut_transcribed.bed
        # And next T>N mutations in genes falling on the negative strand
        ${bedtools_exec} intersect -a ${puva_dir}quintile${q}_expr_genes_skin_negStrand.bed -b ${puva_dir}${sample}_PUVA_mutations_refT.bed >> ${puva_dir}${sample}_Tmut_transcribed.bed

        ## Next find the mutations occurring on the non-transcribed strand. 
        ${bedtools_exec} intersect -a ${puva_dir}quintile${q}_expr_genes_skin_posStrand.bed -b ${puva_dir}${sample}_PUVA_mutations_refT.bed >> ${puva_dir}${sample}_Amut_transcribed.bed

        ${bedtools_exec} intersect -a ${puva_dir}quintile${q}_expr_genes_skin_negStrand.bed -b ${puva_dir}${sample}_PUVA_mutations_refA.bed >> ${puva_dir}${sample}_Amut_transcribed.bed

    done < ${puva_dir}sample_list.txt

    cat ${puva_dir}*_Tmut_transcribed.bed > ${puva_dir}quintile${q}_Tmut_transc.bed
    cat ${puva_dir}*_Amut_transcribed.bed > ${puva_dir}quintile${q}_Amut_transc.bed

done



exit $?