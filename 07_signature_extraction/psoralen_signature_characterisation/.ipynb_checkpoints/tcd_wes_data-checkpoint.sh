#!/bin/bash

## The purpose of this script is to try doing the tcd analysis in puva_tcd.sh and tcd_negative_controls.sh on
## the WES data. The trouble with that may be that the intergenic regions are not covered by the bait set. 
## It turns out that (as expected) this analysis doesn't work. The mutation rates are not comparable in transcribed vs 
## intergenic regions because the intergenic regions are not in the bait set. 

#bsub -o tcd.out -e tcd.err -q normal -R"select[mem>5000] rusage[mem=5000]" -M5000 "bash  ${pso_nfs}07_signature_extraction/strand_asymmetries/tcd_negative_controls.sh"

bedtools_exec=/software/team152/bedtools2/bin/bedtools
puva_dir=/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/
reference=/lustre/scratch117/core/sciops_repository/references/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/all/fasta/Homo_sapiens.GRCh38_full_analysis_set_plus_decoy_hla.fa

    ## Find the mutations occurring on the transcribed strand. 
    ## First G>A and G>T mutations in genes falling on the positive strand
    ${bedtools_exec} intersect -a ${puva_dir}top25pc_expr_genes_skin_posStrand.bed -b ${puva_dir}wes_UV_mutations_refG.bed >> ${puva_dir}wes_Guv_transcribed.bed
    
    # And next C>T and C>A mutations in genes falling on the negative strand
    ${bedtools_exec} intersect -a ${puva_dir}top25pc_expr_genes_skin_negStrand.bed -b ${puva_dir}wes_UV_mutations_refC.bed >> ${puva_dir}wes_Cuv_transcribed.bed

    ## Next find the mutations occurring on the non-transcribed strand.
    ## C>T and C>A mutations on the positive strand
    ${bedtools_exec} intersect -a ${puva_dir}top25pc_expr_genes_skin_posStrand.bed -b ${puva_dir}wes_UV_mutations_refC.bed >> ${puva_dir}wes_Cuv_transcribed.bed

    ${bedtools_exec} intersect -a ${puva_dir}top25pc_expr_genes_skin_negStrand.bed -b ${puva_dir}wes_UV_mutations_refG.bed >> ${puva_dir}wes_Guv_transcribed.bed



exit $?