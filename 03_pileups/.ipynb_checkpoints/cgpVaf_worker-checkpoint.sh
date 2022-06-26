#!/bin/bash


BAM_DIR=$1
reference=$2
output_dir=$3
normal_sample=$4
sample_list=$5
bed_file=$6
mutationType=$7             ## snp/indel

samples=$( echo ${sample_list} | sed 's/,/ /g' )

#declare -a arr=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")
echo "perl-5.16.3 -I /software/CGP/projects/vafCorrect/lib/perl5 /software/CGP/projects/vafCorrect/bin/cgpVaf.pl \
-d ${BAM_DIR} -o ${output_dir} -a ${mutationType} -g ${reference} -be .sample.dupmarked.bam -bo 1 \
-b ${bed_file} -nn ${normal_sample} -tn ${samples} -ct 1 \
-hdr /lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5/shared/ucscHiDepth_0.01_mrg1000_no_exon_coreChrs.bed.gz"


perl-5.16.3 -I /software/CGP/projects/vafCorrect/lib/perl5 /software/CGP/projects/vafCorrect/bin/cgpVaf.pl \
-d ${BAM_DIR} -o ${output_dir} -a ${mutationType} -g ${reference} -be .sample.dupmarked.bam -bo 1 \
-b ${bed_file} -nn ${normal_sample} -tn ${samples} -ct 1 \
-hdr /lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5/shared/ucscHiDepth_0.01_mrg1000_no_exon_coreChrs.bed.gz



exit $?