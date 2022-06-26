#!/bin/bash

sample_list=$1
output_dir=$2
variant_list_dir=$3
verifyBamID_exec=$4
bamDir=$5
REFERENCE=$6


PD=$( awk -v jI="$LSB_JOBINDEX" 'NR==jI {print $1; exit}' < ${sample_list} )
project_nr=$( awk -v jI="$LSB_JOBINDEX" 'NR==jI {print $2; exit}' < ${sample_list} )


${verifyBamID_exec} \
   --Epsilon   1e-12 \
   --Output    ${output_dir}${PD}.verifyBamID.cont \
   --UDPath    ${variant_list_dir}1000g.phase3.10k.b38.exome.vcf.gz.dat.UD  \
   --BamFile   ${bamDir}${project_nr}/${PD}/${PD}.sample.dupmarked.bam \
   --BedPath   ${variant_list_dir}1000g.phase3.10k.b38.exome.vcf.gz.dat.bed  \
   --MeanPath  ${variant_list_dir}1000g.phase3.10k.b38.exome.vcf.gz.dat.mu  \
   --Reference ${REFERENCE}



exit $?