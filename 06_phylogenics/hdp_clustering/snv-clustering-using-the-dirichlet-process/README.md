This set of scripts is intended for running a Dirichlet process to find clusters of SNVs. The clustering output can then be used for downstream analyses including inference of phylogenetic relationships.

There are five main scripts required to run the algorithm, where the final output includes diagnostic plots, cluster assignments of each mutation, and a summary heatmap representation of all SNV clusters found.
Run_Dirichlet_clustering_posthoc.R calls functions from Mutation_Cluster_Dirichlet_Gibbs_sampler_functions.R and Mutation_Cluster_Dirichlet_posthoc_merge_fx.R, which in turn calls a modified version of the cluster assignment correction code stored in label.switching.v2.R.

I have been periodically tweaking these scripts to handle corner cases that the original script did not account for, particularly in cases where clustering is attempted with very small numbers of samples (e.g., ≤ 5). Please check this repository for updates.
Tweaks to the label switching code is speeds up the post processing immensely and is configurable based on the amount of CPU cores and memory available on the machine clustering is performed. I recommend initially setting cores=3 at the top of the script to prevent out of memory errors. The number of execution threads that is most appropriate will depend on the number of mutations and samples analyzed, on top of memory and CPU constraints of the machine.

The general form of the command to start running the algorithm is as follows:
/usr/bin/R < /home/ubuntu/clustering/Run_Dirichlet_clustering_posthoc.R --vanilla --args  \
/directory/containing/clustering/scripts  \
target_alt.csv  \
target_depth.csv  \
mut_contexts.txt  \
/output/directory \
10000 # the number of burn-in iterations to perform before sampling

Below are example input file formats needed to run the algorithm:

target_alt.csv contains alternative allele counts for example from the output of alleleCounter or cgpVaf

# target_alt.csv
chrom,pos,ref,alt,lo006,lo012,lo017,lo042,lo0049,coord_id,mut_id

1,100722726,T,C,0,0,0,7,0,1_100722726,1_100722726_T_C
1,102390132,A,G,0,0,0,7,0,1_102390132,1_102390132_A_G
1,102892042,T,A,18,0,6,2,0,1_102892042,1_102892042_T_A
1,104469757,G,A,6,0,5,0,0,1_104469757,1_104469757_G_A
1,104640755,C,T,7,0,0,0,0,1_104640755,1_104640755_C_T
1,104692843,T,A,0,0,0,13,0,1_104692843,1_104692843_T_A
1,106619682,G,C,11,0,1,4,0,1_106619682,1_106619682_G_C
1,106687053,T,C,0,0,0,8,0,1_106687053,1_106687053_T_C
1,107723952,G,T,14,0,11,2,0,1_107723952,1_107723952_G_T

similarly, target_depth.csv contains depth of coverage information

# target_depth.csv
chrom,pos,ref,alt,coord_id,mut_id,lo006,lo012,lo017,lo042,lo0049

1,100722726,T,C,1_100722726,1_100722726_T_C,85,62,85,103,40
1,102390132,A,G,1_102390132,1_102390132_A_G,85,68,82,85,42
1,102892042,T,A,1_102892042,1_102892042_T_A,105,67,114,86,45
1,104469757,G,A,1_104469757,1_104469757_G_A,78,50,84,77,40
1,104640755,C,T,1_104640755,1_104640755_C_T,100,51,81,75,26
1,104692843,T,A,1_104692843,1_104692843_T_A,94,59,89,91,41
1,106619682,G,C,1_106619682,1_106619682_G_C,97,56,67,64,27
1,106687053,T,C,1_106687053,1_106687053_T_C,99,75,87,83,41
1,107723952,G,T,1_107723952,1_107723952_G_T,94,75,113,79,35

mut_contexts.csv contains trinucleotide context of each variant, which is the output of Context_pull_build37.pl attached (see instructions within script).

# mut_contexts.txt
chrom  pos    ref    alt    WT     CONTEXT GENE   STRAND

1      100722726     T      C      T      AATGCTGCCATGAACATTCAT No_gene No_gene

1      102390132     A      G      A      TGCTTTCTTCATTAATTCTCA ENST00000370103      -1

1      102892042     T      A      T      GAAAAAATGCTCATCATCACT No_gene No_gene

1      104469757     G      A      G      GGATAGGCATGGCACCTGGTG ENST00000441083      1

1      104640755     C      T      C      AGTGCAGCATCAACAATTGTT No_gene No_gene

1      104692843     T      A      T      TATGAATATCTTTCCTTTTCT No_gene No_gene

1      106619682     G      C      G      CTACCATGTAGTGAGAGATGA ENST00000437803      1

1      106687053     T      C      T      CAAATTATCATAGACTGTGTG No_gene No_gene

1      107723952     G      T      G      GCACACCCTTGGTTCACATCT ENST00000370076      1
