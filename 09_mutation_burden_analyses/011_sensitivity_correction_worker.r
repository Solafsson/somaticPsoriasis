
parameters <- commandArgs(TRUE)
patient <- parameters[1]


output_dir="/lustre/scratch119/humgen/projects/psoriasis/burden_analyses/"
cluster_dir="/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/cluster_and_samples/"
binomial_filter_dir="/lustre/scratch119/humgen/projects/psoriasis/binomial_filters/"
hdp_dir="/lustre/scratch119/humgen/projects/psoriasis/phylogenics/hdp_clustering/"

microd_metadata <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/Microdissection_metaData.txt", h=T)
coverage <- microd_metadata[,c("SampleID", "MedianCoverage")]

clusters <- read.csv(paste(cluster_dir, patient, "_cluster_and_samples.csv", sep=""), h=T)
clusters$cluster_id <- gsub("Cl.", paste(patient, "_", sep=""), clusters$cluster_id)

NR <- read.table(paste(binomial_filter_dir, patient, "/", patient, "_NR_pass_sbs.txt", sep=""))
NV <- read.table(paste(binomial_filter_dir, patient, "/", patient, "_NV_pass_sbs.txt", sep=""))

sbs_assign <- read.table(paste(hdp_dir, patient, "/", "ndp_", patient, "_2021_12_09sbs_assignment.txt", sep=""), h=T)



trunc_binom_adjustment <- function(clusters, NR, NV, coverage, hard_limit=4, patient, sbs_assign) {
  cluster_sensitivity <- numeric()
  for(i in 1:nrow(clusters)) {
    clusterID <- clusters$cluster_id[i]
    
    ## Fetch the mutations that are assigned to this cluster
    thisClust <- apply(sbs_assign[sbs_assign$ClusterID==clusterID,c("Chr","Pos","Ref","Alt")], 1, paste, collapse=":")
    thisClust <- gsub(" ", "", thisClust)
    
    ## Calculate the sensitivity of the cluster in each sample
    sensitivities <- as.numeric()
    for(j in 1:ncol(NR)) {
      ## Get the number of reads 
      NR_tmp <- NR[thisClust,j]
      NV_tmp <- NV[thisClust,j]
      NV_tmp <- NV_tmp[NR_tmp>0]
      NR_tmp <- NR_tmp[NR_tmp>0]     # In rare instances, mutations are assigned to a cluster but have zero coverage in the sample. They have coverage in other samples though
      NV_tmp <- NV_tmp[!is.na(NV_tmp)]
      NR_tmp <- NR_tmp[!is.na(NR_tmp)]     
      
      x <- try(binom_mix(NV_tmp,NR_tmp), T)
      if(class(x)=="try-error") {
        ## Model doesn't converge if there are too few reads reporting mutations in this cluster
        ## Corresponds to a sensitivity of 0 (because the clone isn't present in the sample)
        sensitivities <- c(sensitivities,0)
      } else {
        res = binom_mix(NV_tmp,NR_tmp)
        sampleCov <- coverage$MedianCoverage[coverage$SampleID==colnames(NR)[j]]   ## Coverage of the microbiopsy
        sensitivities <- c(sensitivities, mean(unlist(lapply(rpois(n=100000,lambda=sampleCov),function(x) pbinom(q=hard_limit,size=x,p=res$p,lower.tail = F)))))
        
      }
    }
    frac_missing <- 1-sensitivities
    ## The sensitivity is 1- (1-sensitivity in each sample)
    cluster_sensitivity[i] <- 1-prod(frac_missing)
  }
  df <- data.frame(ClusterID=clusters$cluster_id, Mutations=clusters$no.of.mutations.assigned,Sensitivity=cluster_sensitivity,Mutations_adj=clusters$no.of.mutations.assigned/cluster_sensitivity)
  return(df)
}


## Functions originally from Tim Coorens
##################################

## Define the truncated binomial distribution
dbinomtrunc = function(x, size, prob, minx=4) {
  dbinom(x, size, prob) / pbinom(minx-0.1, size, prob, lower.tail=F)
}

estep = function(x,size,p.vector,prop.vector,ncomp, mode){
  ## p.vector = vector of probabilities for the individual components
  ## prop.vector = vector of proportions for the individual components
  ## ncomp = number of components
  p.mat_estep = matrix(0,ncol=ncomp,nrow=length(x))
  for (i in 1:ncomp){
    if(mode=="Truncated") p.mat_estep[,i]=prop.vector[i]*dbinomtrunc(x,size,prob=p.vector[i])
    if(mode=="Full") p.mat_estep[,i]=prop.vector[i]*dbinom(x,size,prob=p.vector[i])
  }
  norm = rowSums(p.mat_estep) ## normalise the probabilities
  p.mat_estep = p.mat_estep/norm
  LL = sum(log(norm)) ## log-likelihood
  
  ## classification of observations to specific components (too crude?)
  which_clust = rep(1,length(x))
  if(ncomp>1){
    which_clust = apply(p.mat_estep, 1, which.max)
  }
  
  list("posterior"=p.mat_estep,
       "LL"=LL,
       "Which_cluster"=which_clust)
}

## Maximisation step
mstep = function(x,size,e.step){
  # estimate proportions
  prop.vector_temp = colMeans(e.step$posterior)
  # estimate probabilities
  p.vector_temp = colSums(x/size*e.step$posterior) / colSums(e.step$posterior)
  
  list("prop"=prop.vector_temp,
       "p"=p.vector_temp)   
}

## EM algorithm
em.algo = function(x,size,prop.vector_inits,p.vector_inits,maxit=5000,tol=1e-6,nclust,binom_mode){
  ## prop.vector_inits =  initial values for the mixture proportions
  ## p.vector_inits =  initial values for the probabilities 
  
  # Initiate EM
  flag = 0
  e.step = estep(x,size,p.vector = p.vector_inits,prop.vector = prop.vector_inits,ncomp=nclust,mode=binom_mode)
  m.step = mstep(x,size,e.step)
  prop_cur = m.step[["prop"]]
  p_cur = m.step[["p"]]
  cur.LL = e.step[["LL"]]
  LL.vector = e.step[["LL"]]
  
  # Iterate between expectation and maximisation steps
  for (i in 2:maxit){
    e.step = estep(x,size,p.vector = p_cur,prop.vector = prop_cur,ncomp=nclust,mode=binom_mode)
    m.step = mstep(x,size,e.step)
    prop_new = m.step[["prop"]]
    p_new = m.step[["p"]]
    
    LL.vector = c(LL.vector,e.step[["LL"]])
    LL.diff = abs((cur.LL - e.step[["LL"]]))
    which_clust = e.step[["Which_cluster"]]
    # Stop iteration if the difference between the current and new log-likelihood is less than a tolerance level
    if(LL.diff < tol){ flag = 1; break}
    
    # Otherwise continue iteration
    prop_cur = prop_new; p_cur = p_new; cur.LL = e.step[["LL"]]
    
  }
  if(!flag) warning("Didn’t converge\n")
  
  BIC = log(length(x))*nclust*2-2*cur.LL
  AIC = 4*nclust-2*cur.LL
  list("LL"=LL.vector,
       "prop"=prop_cur,
       "p"=p_cur,
       "BIC"=BIC,
       "AIC"=AIC,
       "n"=nclust,
       "Which_cluster"=which_clust)
}

binom_mix = function(x,size,nrange=1:5,criterion="BIC",maxit=5000,tol=1e-6, mode="Full"){
  ## Perform the EM algorithm for different numbers of components
  ## Select best fit using the Bayesian Information Criterion (BIC) 
  ## or the Akaike information criterion (AIC)
  i=1
  results = list()
  BIC_vec = c()
  AIC_vec = c()
  
  for (n in nrange){
    ## Initialise EM algorithm with values from kmeans clustering
    init = kmeans(x/size,n)
    prop_init = init$size/length(x)
    p_init = init$centers
    
    results[[i]] = em.algo(x,size,prop.vector_inits = prop_init,p.vector_inits=p_init,nclust=n,maxit,tol,binom_mode=mode)
    BIC_vec = c(BIC_vec,results[[i]]$BIC)
    AIC_vec = c(AIC_vec,results[[i]]$AIC)
    i=i+1
  }
  if (criterion=="BIC"){
    results[[which.min(BIC_vec)]]$BIC_vec=BIC_vec
    return(results[[which.min(BIC_vec)]])
  }
  if (criterion=="AIC"){
    return(results[[which.min(AIC_vec)]])
  }
}
#-------------------------------------------------


results <- trunc_binom_adjustment(clusters, NR, NV, coverage, 4, patient, sbs_assign)
write.table(results, paste(output_dir, patient, "_cluster_sensitivity.txt", sep=""), quote=F, sep="\t", row.names = F)



