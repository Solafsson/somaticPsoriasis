#!/bin/bash

## The purpose of this script is to measure the PUVA-associated mutation rate across gene-expression bins
## to see how and if gene expression affects the mutation rate. 

#bsub -o expr.out -e expr.err -q normal -R"select[mem>5000] rusage[mem=5000]" -M5000 "bash  ${pso_nfs}07_signature_extraction/strand_asymmetries/PUVA_expression.sh"

bedtools_exec=/software/team152/bedtools2/bin/bedtools
puva_dir=/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/
reference=/lustre/scratch117/core/sciops_repository/references/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/all/fasta/Homo_sapiens.GRCh38_full_analysis_set_plus_decoy_hla.fa

for i in {1..10}; do
     awk -v c="$i" '$6=="Bin"c {OFS="\t"; print}' < ${puva_dir}GTex_expression_bins.bed > ${puva_dir}GTex_expression_bin${i}.bed
     rm -f ${puva_dir}All_samples_expr_bin${i}.bed
done


while read sample; do
    for i in {1..10}; do
        ${bedtools_exec} intersect -a ${puva_dir}GTex_expression_bin${i}.bed -b ${puva_dir}${sample}_PUVA_mutations.bed >> ${puva_dir}All_samples_expr_bin${i}.bed
    done
done < ${puva_dir}sample_list.txt

## Count the number of AT/TA dinucleotide sites in each bin:
for i in {1..10}; do
${bedtools_exec} nuc -fi ${reference} -bed ${puva_dir}GTex_expression_bin${i}.bed \
   -pattern AT -C > ${puva_dir}GTex_expression_bin${i}_AT_counts.txt

${bedtools_exec} nuc -fi ${reference} -bed ${puva_dir}GTex_expression_bin${i}.bed \
    -pattern TA -C > ${puva_dir}GTex_expression_bin${i}_TA_counts.txt
done 

## Calculate the number of mutations per AT/TA dinucleotide
rm -f ${puva_dir}Mutation_rate_pr_expr_bin.txt
echo "Bin Mutation_count AT_sites TA_sites Mutation_Rate" > ${puva_dir}Mutation_rate_pr_expr_bin.txt
for i in {1..10}; do
    muts=$( wc -l ${puva_dir}All_samples_expr_bin${i}.bed | awk '{print $1}' )
    AT_sites=$( awk -F'\t' 'BEGIN{SUM=0} NR>1 { SUM=SUM+$NF }END{print SUM}' < ${puva_dir}GTex_expression_bin${i}_AT_counts.txt )
    TA_sites=$( awk -F'\t' 'BEGIN{SUM=0} NR>1 { SUM=SUM+$NF}END{print SUM}' < ${puva_dir}GTex_expression_bin${i}_TA_counts.txt )
    awk -v M="${muts}" -v AT="${AT_sites}" -v TA="${TA_sites}" -v c="$i" 'BEGIN {print c, M, AT, TA, M/(AT+TA)}'  >> ${puva_dir}Mutation_rate_pr_expr_bin.txt
done


exit $?