#!/bin/bash

## The purpose of this script is to prepare the input files for the HDP clustering of mutations. 
#bsub -q normal -o prepare_input_log -e prepare_input_err -R"select[mem>2000] rusage[mem=2000]" -M2000 "bash /nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/06_phylogenics/hdp_clustering/00_prepare_hdp_input.sh"

binomial_dir=/lustre/scratch119/humgen/projects/psoriasis/binomial_filters/
hdp_dir=/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/
script_dir=/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/06_phylogenics/hdp_clustering/

sample_metaData=/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/sample_info/sample_meta.txt
patient_subset=patient79

# Create a list of subjects
awk '( $8=="WGS" || $8=="WES" ) && ($16=="NA" || $16=="Technical_Duplicate") {print $6}' < ${sample_metaData} | sort -u > ${hdp_dir}patientList_hdp.tmp
if [ ! -z "$patient_subset" ]; then
    echo ${patient_subset} | awk -F, -v OFS="\n" '{$1=$1; print}' > ${hdp_dir}patient_subset.tmp
    grep -w -f ${hdp_dir}patient_subset.tmp ${hdp_dir}patientList_binomialFiltering.tmp > ${hdp_dir}tmp
    mv ${hdp_dir}tmp ${hdp_dir}patientList_hdp.tmp
fi

while read patient; do

    mkdir -p ${hdp_dir}${patient}/
    ## Create the input for GRCh38_context_pull.R
    awk ' {print $1}' < ${binomial_dir}${patient}/${patient}_genotype_sbs.txt | awk 'BEGIN{OFS="\t";FS=":"; print "chrom", "pos", "ref", "alt"} NR>1 {FS=":"; OFS="\t"; print $1, $2, $3, $4}' > ${hdp_dir}${patient}/${patient}_context_pull.txt

    ## Create the ndp_alt and ndp_depth files used by Run_Dirichlet_clustering_posthoc.R
    awk 'NR==1 {out="chrom:pos:ref:alt:coord_id:mut_id"; for(i=1; i<=NF; i++) {out=out":"$i} print out} NR>1 {out=$1; for(i=2; i<=NF; i++) {out=out":"$i} print out}' < ${binomial_dir}${patient}/${patient}_NV_pass_sbs.txt | awk 'BEGIN {FS=":"; OFS=":"} NR==1 {print} NR>1 {out=$1":"$2":"$3":"$4":"$1"_"$2":"$1"_"$2"_"$3"_"$4; for(i=5; i<=NF; i++) {out=out":"$i} print out}' | awk 'BEGIN {FS=":"; OFS=","} $1=$1 {print}' > ${hdp_dir}${patient}/${patient}_ndp_alt_bb_flt.csv

    awk 'NR==1 {out="chrom:pos:ref:alt:coord_id:mut_id"; for(i=1; i<=NF; i++) {out=out":"$i} print out} NR>1 {out=$1; for(i=2; i<=NF; i++) {out=out":"$i} print out}' < ${binomial_dir}${patient}/${patient}_NR_pass_sbs.txt | awk 'BEGIN {FS=":"; OFS=":"} NR==1 {print} NR>1 {out=$1":"$2":"$3":"$4":"$1"_"$2":"$1"_"$2"_"$3"_"$4; for(i=5; i<=NF; i++) {out=out":"$i} print out}' | awk 'BEGIN {FS=":"; OFS=","} $1=$1 {print}' > ${hdp_dir}${patient}/${patient}_ndp_depth_bb_flt.csv
    
    /software/R-4.1.0/bin/Rscript ${script_dir}01_GRCh38_context_pull.R ${hdp_dir} ${patient}

done < ${hdp_dir}patientList_hdp.tmp