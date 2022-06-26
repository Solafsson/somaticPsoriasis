#!/bin/bash

## Usage: This script is for making a consensus phylogenic tree for all the sample from an individual. All the steps in
##        the script can be run individually, but this brings it all together for convenience.


## The steps in the tree making process are:
## 1. Read in the binary genotype matrix from beta-binomial filters script.
## 2. Remove samples which are duplicates or should be excluded for any reason. Make the tree with and without.
## 3. Make a fasta file with all the sequences.
## 4. Run MPBoot. Do 3 and 4 separately for subs and indels.
## 5. Assign mutations to the tree branches. Write out the mutations for each branch for signature extraction.
## 6. Adjust the branch lengths

## Post Signature extraction - will be a separate script:
## 7. Plot trees with signatures
## 8. Make trees with branch lengths scaled to adjusted values of signature 1.

## Input: I expect the pileups to have been run and I expect the binomial_dir to contain a number of files...

## for i in 01 03 04 05 06 07 09 11 13 17 23 29 31 32 34 02 10 12 14 16 18 19 20 21 22 24 26 27 28 30 08 115 37 46 55 76 97 38 47 56 79 98 109 29 39 48 57 99 102 40 49 58 80 103 110 41 50 60 84 104 42 51 62 86 112 25 43 52 70 87 105 106 113 36 44 53 95 111; do bsub -o logs/patient${i}_snv.out -e logs/patient${i}_snv.err -M1000 -R'span[hosts=1] select[mem>1000] rusage[mem=1000]' "bash /nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/06_phylogenics/make_MPBoot_tree_master.sh patient${i} snv"; done

# bsub -o logs/patient34_WGS_snv.out -e logs/patient34_WGS_snv.err -M1000 -R'span[hosts=1] select[mem>1000] rusage[mem=1000]' "bash /nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/06_phylogenics/make_MPBoot_tree_master.sh patient34_WGS snv"

pileup_dir=/lustre/scratch119/humgen/projects/psoriasis/pileups/
binomial_dir=/lustre/scratch119/humgen/projects/psoriasis/binomial_filters/
script_dir=/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/06_phylogenics/

patient=$1
mutType=$2        ## snv or indel
MPBoot_dir=/lustre/scratch119/humgen/projects/psoriasis/phylogenics/MPBoot/${mutType}/


mkdir -p ${MPBoot_dir}consensus_trees
mkdir -p ${MPBoot_dir}input_fasta_files
mkdir -p ${MPBoot_dir}binary_gt_matrices
mkdir -p ${MPBoot_dir}branch_mut_assignment
mkdir -p ${MPBoot_dir}tree_plots
mkdir -p ${MPBoot_dir}adj_mutCounts


if [ "$mutType" == "snv" ]; then
    echo "Making fasta file:"
    echo "/software/R-3.6.1/bin/Rscript ${script_dir}make_fasta_and_binary_gt_matrix.r ${patient} ${binomial_dir}${patient}/${patient}_genotype_allSubs.txt ${MPBoot_dir}input_fasta_files/"
    
    ## Currently only using the SBS mutations because the sequences must all have the same length and a different number of
    ## DBS mutations is messing up the script. Think about how to incorporate them and use _genotype_allSubs.txt
    /software/R-3.6.1/bin/Rscript ${script_dir}make_fasta_and_binary_gt_matrix.r ${patient} ${mutType} \
    ${binomial_dir}${patient}/${patient}_genotype_sbs.txt ${MPBoot_dir}input_fasta_files/

    echo "Running MPBoot"
    echo "${script_dir}run_MPBoot.sh ${MPBoot_dir} ${MPBoot_dir}input_fasta_files/${patient}_trees.fa ${patient}"
    ${script_dir}run_MPBoot.sh ${MPBoot_dir} ${MPBoot_dir}consensus_trees/ ${MPBoot_dir}input_fasta_files/${patient}_trees.fa ${patient}

    echo "Processing trees:"
    echo "/software/R-3.6.1/bin/Rscript ${script_dir}process_trees.R ${patient} ${binomial_dir}${patient}/${patient}_genotype_allSubs.txt ${MPBoot_dir}consensus_trees/${patient}/${patient}.treefile ${binomial_dir}${patient}/${patient}_NR_pass_allSubs.txt ${binomial_dir}${patient}/${patient}_NV_pass_allSubs.txt ${binomial_dir}${patient}/${patient}_medVAF_and_Unadj_MutCount.txt ${MPBoot_dir}"

    /software/R-4.1.0/bin/Rscript ${script_dir}process_trees.r ${patient} ${binomial_dir}${patient}/${patient}_genotype_allSubs.txt \
    ${MPBoot_dir}consensus_trees/${patient}/${patient}.treefile ${binomial_dir}${patient}/${patient}_NR_pass_allSubs.txt \
    ${binomial_dir}${patient}/${patient}_NV_pass_allSubs.txt \
    ${MPBoot_dir} ${mutType}
fi

if [ "$mutType" == "indel" ]; then
    echo "Processing indels for patient ${patient}"

    echo "Making fasta file:"
    /software/R-3.6.1/bin/Rscript ${script_dir}make_fasta_and_binary_gt_matrix.r ${patient} ${mutType} \
    ${binomial_dir}${patient}/${patient}_genotype_indels.txt ${MPBoot_dir}input_fasta_files/

    echo "Running MPBoot"
    echo "${script_dir}run_MPBoot.sh ${MPBoot_dir} ${MPBoot_dir}input_fasta_files/${patient}_trees.fa ${patient}"
    ${script_dir}run_MPBoot.sh ${MPBoot_dir} ${MPBoot_dir}consensus_trees/ ${MPBoot_dir}input_fasta_files/${patient}_trees.fa ${patient}

    ## I purposely pass the snv tree to the process_trees.r script as I want to only retain those indels that adhere
    ## to the more reliable snv tree structure.
    echo "Processing trees:"
    /software/R-4.1.0/bin/Rscript ${script_dir}process_trees.r ${patient} ${binomial_dir}${patient}/${patient}_genotype_indels.txt \
    ${MPBoot_dir}consensus_trees/${patient}/${patient}.treefile \
    ${binomial_dir}${patient}/${patient}_NR_pass_indel.txt \
    ${binomial_dir}${patient}/${patient}_NV_pass_indel.txt \
    ${MPBoot_dir} ${mutType}
fi

exit $?