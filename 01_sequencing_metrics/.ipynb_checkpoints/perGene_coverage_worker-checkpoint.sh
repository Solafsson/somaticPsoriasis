#!/bin/bash


sample=$1
projectNr=$2
output_dir=$3
gene_coordinates=$4
gene=$5

data_dir=/nfs/cancer_ref01/nst_links/live/


/software/team152/samtools-1.11/samtools depth -a -b ${gene_coordinates} ${data_dir}${projectNr}/${sample}/${sample}.sample.dupmarked.bam -q 30 -Q 30 -s > ${output_dir}${sample}_${gene}_samtools_depth.txt

awk '{print $3}' < ${output_dir}${sample}_${gene}_samtools_depth.txt | sort -n | awk -v SAM="$sample" 'BEGIN {isZero=0; minThree=0; minFour=0; minFive=0; minTen=0} {if($1==0) {isZero=isZero+1}; if($1<4) {minThree=minThree+1}; if($1<5) {minFour=minFour+1}; if($1<6) {minFive=minFive+1}; if($1<11) {minTen=minTen+1} { a[i++]=$1; }} END { x=int((i+1)/2); if (x < (i+1)/2) print SAM,isZero, minThree, minFour, minFive, minTen, (a[x-1]+a[x])/2; else print SAM, isZero, minThree, minFour, minFive, minTen, a[x-1]; }' >> ${output_dir}samtools_${gene}_coverage.txt

rm ${output_dir}${sample}_${gene}_samtools_depth.txt

done < ${gene_coordinates}



exit $?


