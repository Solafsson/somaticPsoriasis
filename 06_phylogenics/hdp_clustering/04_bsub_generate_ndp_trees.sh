#!/bin/bash

SCRIPT_DIR=/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/06_phylogenics/hdp_clustering/

hdp_dir=/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/ # Dir that contains output of ndp clustering
tree_out_dir=/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/trees/
min_vaf_threshold=0.05 #default 0.10
min_mut_count=10 # default 50


sample_metaData=/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/sample_info/sample_meta.txt
patient_subset=$1

# Create a list of subjects
awk '( $8=="WGS" || $8=="WES" ) && ($16=="NA" || $16=="Technical_Duplicate") {print $6}' < ${sample_metaData} | sort -u > ${hdp_dir}patientList_hdp.tmp
if [ ! -z "$patient_subset" ]; then
    echo ${patient_subset} | awk -F, -v OFS="\n" '{$1=$1; print}' > ${hdp_dir}patient_subset.tmp
    grep -w -f ${hdp_dir}patient_subset.tmp ${hdp_dir}patientList_hdp.tmp > ${hdp_dir}tmp
    mv ${hdp_dir}tmp ${hdp_dir}patientList_hdp.tmp
fi


MEM=4000
CORE=1
mkdir -p $hdp_dir/logfiles
# Failed: 7,110,54,55,79,37
while read patient; do

#cp ${hdp_dir}${patient}/${patient}_ndp_alt_bb_flt.csv  ${hdp_dir}${patient}/ndp_${patient}_2021_12_09/${patient}_bb_pass_snvs_all.csv
cp ${hdp_dir}${patient}/${patient}_ndp_alt_bb_flt.csv  ${hdp_dir}${patient}/ndp_${patient}_2022_02_15/${patient}_bb_pass_snvs_all.csv
# Run R code
  bsub -J ${patient}.ndp_tree -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -n $CORE -q small \
  -e $hdp_dir/logfiles/${patient}.stderr -o $hdp_dir/logfiles/${patient}.stdout \
  "/software/R-4.1.0/bin/R < $SCRIPT_DIR/05_ndp_tree_generation.R --vanilla --args \
          ${patient} \
          ${hdp_dir}${patient}/ndp_${patient}_2021_12_09/ \
          $SCRIPT_DIR \
          $tree_out_dir \
          $min_vaf_threshold \
          $min_mut_count"

done < ${hdp_dir}patientList_hdp.tmp