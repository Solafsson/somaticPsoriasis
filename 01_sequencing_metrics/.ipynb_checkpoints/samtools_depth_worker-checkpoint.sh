#!/bin/bash

## The purpose of this script is to calculate the true coverage of the WES samples as the
## CanApps pipeline gives clearly inflated estimates


sample=$1
projectNr=$2
output_dir=$3

data_dir=/nfs/cancer_ref01/nst_links/live/
baitset_bed=/lustre/scratch119/humgen/projects/psoriasis/resources/Human_all_exon_V5.bed


## See documentation of genomecov here: https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html
#/software/team152/bedtools genomecov -bga -ibam ${data_dir}${projectNr}/${sample}/${sample}.sample.dupmarked.bam > ${output_dir}${sample}.genomecov.txt

#bsub -o test.out -e test.err -R"select[mem>2000] rusage[mem=2000]" -M2000 "/software/team152/bedtools coverage -a ../resources/Human_all_exon_V5.bed -b /nfs/cancer_ref01/nst_links/live/2545/P01H_LL_6/P01H_LL_6.sample.dupmarked.bam -hist > P01H_LL_6.HistCoverageBED.bedgraph"


/software/team152/samtools-1.11/samtools depth -a -b ${baitset_bed} ${data_dir}${projectNr}/${sample}/${sample}.sample.dupmarked.bam -q 30 -Q 30 -s > ${output_dir}${sample}.samtools_depth.txt

## The output file is quite large and I don't want to keep it. Instead, I'm going to extract the number of bases with 0 coverage,
## the number with at least 3, 4,5 and 10X coverage and the median for the sample. Then I'll delete the file.
awk '{print $3}' < ${output_dir}${sample}.samtools_depth.txt | sort -n | awk -v SAM="$sample" 'BEGIN {isZero=0; minThree=0; minFour=0; minFive=0; minTen=0} {if($1==0) {isZero=isZero+1}; if($1<4) {minThree=minThree+1}; if($1<5) {minFour=minFour+1}; if($1<6) {minFive=minFive+1}; if($1<11) {minTen=minTen+1} { a[i++]=$1; }} END { x=int((i+1)/2); if (x < (i+1)/2) print SAM,isZero, minThree, minFour, minFive, minTen, (a[x-1]+a[x])/2; else print SAM, isZero, minThree, minFour, minFive, minTen, a[x-1]; }' >> ${output_dir}samtools_depth_coverage.txt

rm ${output_dir}${sample}.samtools_depth.txt

exit $?

