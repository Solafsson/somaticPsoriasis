#!/bin/bash

## The purpose of this script is to submit the sensitivity correction to the cluster.
## See 011_sensitivity_correction_worker.r
## /nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/09_mutation_burden_analyses/01_sensitivity_correction_master.sh
patient_list=$1

script_dir=/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/09_mutation_burden_analyses/

while read patient; do

    bsub -q normal -o ${patient}_sensitivity_log -e ${patient}_sensitivity_err -R"select[mem>1000] rusage[mem=1000]" -M1000 \
    "/software/R-4.1.0/bin/Rscript ${script_dir}011_sensitivity_correction_worker.r ${patient}"

done < ${patient_list}

exit $?