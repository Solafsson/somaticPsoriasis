## The purpose of this script is to perform HDP-signature extraction for single base substitutions

#.libPaths("/lustre/scratch114/projects/crohns/somatic_ibd_p1/R_packages_farm5_R3.6_install/")
.libPaths("/lustre/scratch119/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")
library(hdp)

parameters <- commandArgs(TRUE)
mut_type <- parameters[1]
output_dir <- parameters[2]
jobIndex <- as.numeric(parameters[3])

if(mut_type=="snv") {
    all_sigs <- read.csv("/lustre/scratch119/humgen/projects/psoriasis/resources/sigProfiler_SBS_signatures_summing_to_one.csv", h=T,row.names = 1, stringsAsFactors = F)
        
    ## Filter the signatures.
# Set prior to those signatures found to be active in most normal tissues, UV-related signatures, signatures found in skin cancers and signatures
# who's activity was acelerated in IBD. 
gdsigs <- c("SBS1", "SBS5","SBS2", "SBS13","SBS7a","SBS7b","SBS7c","SBS7d", "SBS17a", "SBS17b","SBS18","SBS38")
    
} 


load(paste(output_dir, "branch_mutational_matrix.txt", sep=""))

randseed <- 42*(jobIndex) + 1
print(paste0("Job index:", jobIndex))
print(paste0("Random seed:", randseed))


# Format the mutation matrix
tnt <- mut_mat_t
newColNames <- paste(substr(colnames(mut_mat_t), 3,3), ".", substr(colnames(mut_mat_t), 5,5), ".in.", substr(colnames(mut_mat_t), 1,1), substr(colnames(mut_mat_t), 3,3), substr(colnames(mut_mat_t), 7,7), sep="")
colnames(tnt) <- newColNames
patients <- sapply(strsplit(rownames(tnt), "_"), "[[", 1)

# Retain only the specified prior signatures
sigs <- as.matrix(all_sigs[,gdsigs])

## Set up the priors for HDP
prior_pseudoc <- rep(10000, ncol(sigs)) # this is meant to be a vector of pseudocounts contributed by each prior distribution.
test_prior <- hdp_prior_init(sigs, prior_pseudoc, hh=rep(1, 96), alphaa=c(1, 1), alphab=c(1, 1))

# add three more CPs for the data that we will add: 
# one for the branch from the grandparent to the MRCA of the patients, 
# one for all the branches from the MRCA of the patients to each patient
# one to go from each patient to all their child nodes (same cp for all patients)
test_prior <- hdp_addconparam(test_prior, alphaa=c(1,1,1), alphab=c(1,1,1))



# now add the hierarchical structure of the dataset.
# I want: 
# ppindex 0 for the top node. - ALREADY PRESENT
# ppindex 1 for each of the 28 signatures. - ALREADY PRESENT
# ppindex 1 for the MRCA of all the patients. - TO ADD
# ppindex 30 (i.e. the MRCA of all the patients) for all 42 patients. - TO ADD
# and then each branch has the ppindex of the patient that it comes from. - TO ADD

ptsmrcappindex <- 1
ptsppindices <- rep(c(1 + ncol(sigs) + 1), length(unique(patients)))

branchppindices <- as.numeric(unlist(sapply(unique(patients), function(patient) {
  tnum <- length(which(patients==patient))
  tindex <- (1 + ncol(sigs) + 1) + which(unique(patients)==patient)
  rep(tindex, tnum)
})))

newppindices <- c(ptsmrcappindex, ptsppindices, branchppindices)

# For concentration parameters, give one for every level of the hierarchy:
# cpindex 1 for the grandparent. - ALREADY PRESENT
# cpindex 2 for all the signatures - ALREADY PRESENT
# cpindex 3 for the parent of the patients - NEED TO ADD
# cpindex 4 for all the patients - NEED TO ADD
# cpindex 5 for all the patients to their branches - NEED TO ADD
newcpindices <- c(3, rep(4, length(unique(patients))), rep(5, nrow(tnt)))

# add dp nodes: 
# one as the parent of all the patients
# one for each of the 42 patients
# one for every single branch.
# i.e. this should be the same as the number of new ppindicies and cpindices
test_prior <- hdp_adddp(test_prior,
                        numdp = length(newppindices),
                        ppindex = newppindices,
                        cpindex = newcpindices)

# need to make sure that I am adding data in a way that matches the data to the terminal nodes.
# dpindices are the indices of the terminal nodes. The first terminal node is number 73.
test_prior <- hdp_setdata(hdp=test_prior, dpindex=(max(newppindices) + 1:nrow(tnt)), data=tnt)


# run chain
test_pr <- dp_activate(test_prior, dpindex=(1+ncol(sigs)+1):numdp(test_prior), initcc=(ncol(sigs)+round(ncol(sigs)/10)))
test_chlist <- hdp_posterior(test_pr, burnin=100000, n=100, space=2000, cpiter=3, seed=randseed)

# save data.
assign(paste0("hdp_", jobIndex), test_chlist)

save(list=paste0("hdp_", jobIndex), file=paste(output_dir, "hdpout_", jobIndex, ".RData", sep=""))


