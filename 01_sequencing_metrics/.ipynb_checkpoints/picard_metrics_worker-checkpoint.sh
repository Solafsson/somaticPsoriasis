#!/bin/bash

#bsub -o picard.out -e picard.err -R"select[mem>6000] rusage[mem=6000]" -M6000 "/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/05_early_exploratory_analyses/picard_metrics.sh"

sample=$1
projectNr=$2
output_dir=$3

#First, you need to create a list of intervals for the exome. 
# You only need to do this once.
#java -jar /lustre/scratch119/casm/team268im/fa8/bin/picard-tools-1.131/picard.jar BedToIntervalList \
#I=/lustre/scratch119/humgen/projects/psoriasis/resources/Human_all_exon_V5.bed \
#O=/lustre/scratch119/humgen/projects/psoriasis/coverage_calculation/Human_exome_v5_Covered.intervals \
#SD=/lustre/scratch117/core/sciops_repository/references/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/all/picard/Homo_sapiens.GRCh38_full_analysis_set_plus_decoy_hla.fa.dict


java -jar -Xmx5g /lustre/scratch119/casm/team268im/fa8/bin/picard-tools-1.131/picard.jar CalculateHsMetrics \
INPUT=/nfs/cancer_ref01/nst_links/live/${projectNr}/${sample}/${sample}.sample.dupmarked.bam \
OUTPUT=${output_dir}${sample}.all.hsmetrics \
TARGET_INTERVALS=/lustre/scratch119/humgen/projects/psoriasis/coverage_calculation/Human_exome_v5_Covered.intervals \
BAIT_INTERVALS=/lustre/scratch119/humgen/projects/psoriasis/coverage_calculation/Human_exome_v5_Covered.intervals \
VALIDATION_STRINGENCY=SILENT 



# echo "SampleID MEAN_TARGET_COVERAGE PCT_PF_UQ_READS PCT_OFF_BAIT" > coverage_report.txt
# while read sample; do awk -v SAM="$sample" '$1=="Human_exome_v5_Covered" {print SAM, $22, $10,$19}' < ${sample}.all.hsmetrics >> coverage_report.txt; done < <( awk '$16=="NA" && $8=="WES" {print $1}' < /nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/sample_info/sample_meta.txt )



exit $?