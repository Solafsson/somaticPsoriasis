#!/bin/bash

## The purpose of this script is to submit the R-script running the pentanucleotide model 
## for pathways to the cluster

working_dir=/lustre/scratch119/humgen/projects/psoriasis/selection_analyses/pathway_pentamodel/

while read geneList; do
bsub -o ${working_dir}pentanucleotide_${geneList}.out -e ${working_dir}pentanucleotide_${geneList}.err -q long -R"select[mem>20000] rusage[mem=20000]" -M20000 "/software/R-4.1.0/bin/Rscript /nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/08_selection_analyses/pathway_dnds_pentanuc_model.r ${geneList}" 

done < /lustre/scratch119/humgen/projects/psoriasis/selection_analyses/pathway_dNdS_geneLists.txt

exit $?