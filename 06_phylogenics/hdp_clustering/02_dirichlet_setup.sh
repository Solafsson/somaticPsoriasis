#!/bin/bash

## Example:
#bsub -q normal -o dirichlet_setup_log -e dirichlet_setup_err -R"select[mem>2000] rusage[mem=2000]" -M2000 "bash /nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/06_phylogenics/hdp_clustering/02_dirichlet_setup.sh"


script_dir=/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/06_phylogenics/hdp_clustering/
hdp_dir=/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/

sample_metaData=/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/sample_info/sample_meta.txt
patient_subset=$1

# Create a list of subjects
awk '( $8=="WGS" || $8=="WES" ) && ($16=="NA" || $16=="Technical_Duplicate") {print $6}' < ${sample_metaData} | sort -u > ${hdp_dir}patientList_hdp.tmp
if [ ! -z "$patient_subset" ]; then
    echo ${patient_subset} | awk -F, -v OFS="\n" '{$1=$1; print}' > ${hdp_dir}patient_subset.tmp
    grep -w -f ${hdp_dir}patient_subset.tmp ${hdp_dir}patientList_hdp.tmp > ${hdp_dir}tmp
    mv ${hdp_dir}tmp ${hdp_dir}patientList_hdp.tmp
fi

# Change to working dir (where error logs will go)
mkdir -p $hdp_dir/logfiles

#source ~/.Renviron


burnin=15000 # reccomend 15,000
MEM=6000
CORE=8
q=basement

while read patient; do

# Run R code
  bsub -J dirichlet_setup_${patient} -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -n $CORE -q $q \
          -e $hdp_dir/logfiles/$patient.stderr -o $hdp_dir/logfiles/$patient.stdout \
          "/software/R-4.1.0/bin/R < ${script_dir}snv-clustering-using-the-dirichlet-process/Run_Dirichlet_clustering_posthoc.R  --vanilla --args ${script_dir}snv-clustering-using-the-dirichlet-process/ ${hdp_dir}${patient}/$patient*ndp_alt_bb_flt.csv ${hdp_dir}${patient}/$patient*ndp_depth_bb_flt.csv ${hdp_dir}${patient}/$patient*mut_context_GRCh38.txt  ${hdp_dir}${patient}/ $burnin"
echo ${patient}

done < ${hdp_dir}patientList_hdp.tmp
