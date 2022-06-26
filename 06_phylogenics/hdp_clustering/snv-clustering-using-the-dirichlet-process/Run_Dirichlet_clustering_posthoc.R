# Run Dirichlet clustering
# --
# Runs Dirichlet clustering and label switching on low input mutation data
#
# Args (command line):
#   1: Path to git repository (ie. parent directory where repository was cloned into)
#   2: Path to CSV containing counts of alt base (rows are variants, any columns starting with 'PD' are samples)
#   3: Path to CSV containing depth (rows are variants, any columns starting with 'PD' are samples)
#   4: Path to table with mutation contexts (generated with Peter Campbell's "Context_pull_build37.pl" script)
#   5: Directory to output results
#   6: Number of burnin cycles, eg 10000
# 
# If command line args are not provided, default file locations are used.
# 
# --
# /// Author --- PETER CAMPBELL
# /// Edits --- Simon Brunner
# /// Edit date --- 11-MAY-2018
# /// Edits --- Stanley Ng
# /// Edit date --- 25-MAY-2019

#############################################
## -- Receive command line args --
args = commandArgs(trailingOnly=TRUE)
args

# Check for command line args, else set defaults
if (length(args)==0) {
  repo_location = '/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/06_phylogenics/hdp_clustering/snv-clustering-using-the-dirichlet-process/'
  fpath_alt_csv = '/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/patient07_ndp_alt_bb_flt.csv'
  fpath_depth_csv = '/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/patient07_ndp_depth_bb_flt.csv'
  fpath_mutant_locations = '/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/patient07_mut_context_GRCh38.txt'
  output_dir = '/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/'
  num_burnin = 10000
} else if (length(args)==6) {
  repo_location = args[1]
  fpath_alt_csv = args[2]
  fpath_depth_csv = args[3]
  fpath_mutant_locations = args[4]
  output_dir = args[5]
  num_burnin = as.numeric(args[6])
} else {
  stop("Script takes either exactly 0 or 6 arguments.", call.=FALSE)
}



# Report arguments
message(sprintf('Reading alt counts from %s', fpath_alt_csv))
message(sprintf('Reading depth from %s', fpath_depth_csv))
message(sprintf('Reading mutation contexts from %s', fpath_mutant_locations))
message(sprintf('Writing output to %s', output_dir))
message(sprintf('Number of burnin cycles is %s', num_burnin))

#############################################
## -- Library loading and sourcing --

.libPaths("/lustre/scratch119/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")

library(RColorBrewer)
library(label.switching)
library(philentropy)
library(ggplot2)
library(GGally)
library(dplyr)
library(data.table)
library(here)
library(farver)

# source(file.path(repo_location, 'lib', 'Mutation_Cluster_Dirichlet_Gibbs_sampler_functions.R'))
# source(file.path(repo_location, 'lib', 'Mutation_Cluster_Dirichlet_posthoc_merge_fx.R'))
source(file.path(repo_location,'Mutation_Cluster_Dirichlet_Gibbs_sampler_functions.R'))
source(file.path(repo_location,'Mutation_Cluster_Dirichlet_posthoc_merge_fx.R'))

today = format(Sys.Date(), "%Y_%m_%d")
#this_hash = paste(sample(letters,5,TRUE),collapse="")
output_dir = file.path(output_dir, paste('ndp', unlist(strsplit(basename(fpath_alt_csv), split="_"))[1], today,sep="_"))

dir.create(output_dir, showWarnings = F)

# Define params
number.iterations <- 25000
number.burn.in <- 15000
number.burn.in <- num_burnin

set.seed(28)

#############################################
## -- Run clustering --

mut.spec <- read.table(fpath_mutant_locations, sep="\t", header=TRUE, 
                       stringsAsFactors = FALSE)
mut.spec <- mut.spec[ mut.spec$chrom != "chrY",]

lymph.var <- read.table(fpath_alt_csv, sep=",", header=TRUE, stringsAsFactors = FALSE)
lymph.var <- lymph.var[lymph.var$chrom != "chrY",] # Remove Y chrom muts
lymph.y <- as.matrix(lymph.var[,which(grepl('P', names(lymph.var)))])

lymph.depth <- read.table(fpath_depth_csv, sep=",", header=TRUE, stringsAsFactors = FALSE)
lymph.depth <- lymph.depth[lymph.depth$chrom != "chrY",]
lymph.N <- as.matrix(lymph.depth[,which(grepl('P', names(lymph.depth)))])

lymph.gs <- muts.cluster.dirichlet.gibbs(C = 100, y = lymph.y, N = lymph.N, iter=number.iterations,
                                         iter.chunk=100, outdir=output_dir)


short.names <- names(lymph.var)[which(grepl('P', names(lymph.var)))]
#short.names[2:length(short.names)] <- substring(short.names[2:length(short.names)], 8)

save.image(file.path(output_dir, 'Rsession.dat'))

#############################################

# Generate convergence plots

pdf(file.path(output_dir, "Convergence_plots.pdf"), paper="a4", width = 8, height = 10)

for (i in c(10, 50, 100, 200, 500, seq(1000, number.iterations, 1000))) {
  heatmap.iter(lymph.gs, curr.iter = i, samp.names = short.names, min.threshold = 10)
}

plot(1:number.iterations, lymph.gs[["alpha"]], type="l", xlab="Iteration", ylab = "Alpha", main = )

par(mfrow = c(4,1))
for (i in 1:4) {
  label.switch.plot(lymph.gs, num.samples = 10)
}

dev.off()

#############################################

# Correct label switching phenomenon

lymph.ls <- correct.label.switch(lymph.gs, method.set = c("ECR"), burn.in = number.burn.in)
ls_tbl = tibble(clust_assign = as.numeric(lymph.ls$clusters)) %>%
  mutate(mut_id = row_number())
fwrite(ls_tbl, file.path(output_dir, 'clust_assign.csv'))

saveRDS(lymph.ls, file.path(output_dir, 'object_lymph.ls.dat'))
save.image(file.path(output_dir, 'Rsession_ls.dat'))
