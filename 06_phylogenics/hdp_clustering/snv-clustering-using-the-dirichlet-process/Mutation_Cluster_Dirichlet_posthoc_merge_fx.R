#############################################

# Merge clusters where differences are only related to sequencing errors

merge.clusters.differ.by.seq.error <- function(gs.out, ls.out, min.true.VAF = 0.025, overlap.frac=0.025, min.num.muts = 20, ls.type="ECR", burn.in, version=2) {
  # Function to merge clusters that are only separated from one another by seq error rates
  # Need for the function is that some clusters are falsely split by the level of seq error
  ## in samples where the mutations are absent (ie, some base positions noisier than others)
  
  # Two versions of this algorithm:
  # Version 1 = For 2 clusters, looks to see whether iterations of MCMC alternate sign for each non-zero sample
  # Version 2 = For 2 clusters, looks to see whether mean locations (pi.h.j) for each non-zero sample differ by less that a certain threshold 
  
  # gs.out is the output from the Gibbs sampler
  # ls.out is the output from the label.switching() function
  # min.true.VAF is the VAF threshold at which a sample will be counted as carrying the mutations
  # overlap.frac (VERSION 1) is the minimum fraction of MCMC estimates of pi that overlap in the two clusters for each non-zero sample required for merging - increase its value for less merging; decrease its value for more merging 
  # overlap.frac (VERSION 2) is the maximum difference between pi.h.j estimates for each non-zero sample
  # min.num.muts is the minimum number of mutations in the clusters to worry about merging them
  
  ls.final <- ls.out$clusters[1,]
  unique.clusters <- table(ls.final)
  kept_clusters = which(unique.clusters >= min.num.muts)
  unique.clusters <- unique.clusters[kept_clusters]
  unique.clusters <- names(unique.clusters)[order(unique.clusters, decreasing = TRUE)]
  
  y <- gs.out[["y1"]]
  N <- gs.out[["N1"]]
  merge.spl <- gs.out[["merge.split"]]
  
  ls.perm <- ls.out$permutations[ls.type][[1]]
  pi.h.j <- gs.out[["pi.h.j"]]
  num.clusters <- dim(pi.h.j)[3]
  num.samples <- dim(pi.h.j)[2]
  num.iter <- dim(pi.h.j)[1]
  
  pi.mcmc <- array(NA, dim = c(num.iter - burn.in, num.clusters, num.samples))
  for (i in (burn.in+1):num.iter) {pi.mcmc[i-burn.in,,] <- t(pi.h.j[i,,])}
  # Permute according to label-switching() defined optimal labelling 
  pi.mcmc <- permute.mcmc(mcmc = pi.mcmc, permutations = ls.perm)$output
  
  # For successful merge-split steps after burn.in, turn estimates of pi before the merge-split to NA
  merge.spl <- merge.spl[(burn.in+1):num.iter,]
  merge.spl <- merge.spl[merge.spl$decision,]

  # condition to handle cases where no merging of clusters is required
  if(dim(merge.spl)[1] > 0){

  merge.spl$new.cluster.loc <- as.numeric(merge.spl$new.cluster.loc)
  for (i in 1:nrow(merge.spl)) {
    final.cluster.1 <- ls.perm[merge.spl$iter[i] - burn.in, merge.spl$first.cluster[i]]
    final.cluster.2 <- ls.perm[merge.spl$iter[i] - burn.in, merge.spl$second.cluster[i]]
    final.cluster.3 <- ls.perm[merge.spl$iter[i] - burn.in, merge.spl$new.cluster.loc[i]]
    
    # instead of assigning N/A's, assign 0's instead to ensure no N/A's are present
    # during subsequent mean.vaf calculation and heatmap plotting ???
    if (merge.spl$first.cluster[i] == merge.spl$second.cluster[i]) {
      pi.mcmc[1:(merge.spl$iter[i] - burn.in), final.cluster.1, ] <- NA
      pi.mcmc[1:(merge.spl$iter[i] - burn.in), final.cluster.3, ] <- NA
    } else {
      pi.mcmc[1:(merge.spl$iter[i] - burn.in), final.cluster.1, ] <- NA
      pi.mcmc[1:(merge.spl$iter[i] - burn.in), final.cluster.2, ] <- NA
    }
  }

  }
  
  # since above can assign some N/A's to pi.mcmc some columns in mean.vaf will be N/A,
  # should N/A's be converted to zeros, or should merge step be skipped if N/A's detected ???
  mean.vaf <- apply(pi.mcmc, MARGIN = 2, FUN = colMeans, na.rm=TRUE)
  if (version == 2) {
    S.i <- read.table(gs.out[["S.i.filename"]], header=TRUE, sep="\t", stringsAsFactors = FALSE,
                      skip = burn.in)
    size.clusters <- t(apply(S.i, MARGIN = 1, 
                             FUN = function(a) {table(factor(a, levels = 1:num.clusters))}))
  }
  
  # Check whether to merge clusters
  # Merge if in all samples with VAFs above the seq error threshold, MCMC pi overlaps > overlap.frac times 
  for (i in 1:(length(unique.clusters)-1)) {
    ind.i <- as.double(unique.clusters[i])

    # to handle the case of N/A values in mean.vaf columns:
    # skip merge step if N/A values detected
    if(length(is.na(mean.vaf[,ind.i])) > 0){next}

    for (j in (i+1):length(unique.clusters)) {
      ind.j <- as.double(unique.clusters[j])

      # to handle the case of N/A values in mean.vaf columns:
      # skip merge step if N/A values detected
      if(length(is.na(mean.vaf[,ind.j])) > 0){next}
      # zero out NA's
      # if(length(is.na(mean.vaf[,ind.i])) > 0){
      #  idx = which(is.na(mean.vaf[,ind.i]))
      #  mean.vaf[idx,ind.i] = 0
      # }

      nonzero.samps <- which(mean.vaf[,ind.i] >= min.true.VAF | mean.vaf[,ind.j] >= min.true.VAF )
      num.overlaps <- sapply(1:dim(pi.mcmc)[3], function(a,k,x,y) {
        temp <- na.exclude(a[,x,k] > a[,y,k]);
        if (length(temp)>0) {sum(temp) / length(temp)} else {0}}, 
        a=pi.mcmc, x=ind.i, y=ind.j)
      if (version == 1 & 
          all(num.overlaps[nonzero.samps] > overlap.frac & 
              num.overlaps[nonzero.samps] < 1-overlap.frac)) {
        # Clusters should be merged
        # Merge them into the larger cluster (ind.i)
        ls.final[ls.final == ind.j] <- ind.i
      } else if (version == 2 & all(abs(mean.vaf[,ind.i] - mean.vaf[,ind.j]) < overlap.frac)) {
        # Clusters should be merged - weight the pi estimates by size of cluster
        ls.final[ls.final == ind.j] <- ind.i
        pi.mcmc[ , ind.i, ] <- (size.clusters[, ind.i] * pi.mcmc[ , ind.i, ] + 
                                  size.clusters[, ind.j] * pi.mcmc[ , ind.j, ]) /
          (size.clusters[, ind.i] + size.clusters[, ind.j])
      }
    }
  }
  return(list(final.cluster = ls.final, pi.mcmc = pi.mcmc))
}


#############################################

# Mutational spectra per cluster

context.extract_posthoc <- function(muts, ls.merge.out, plot.heatmap = TRUE, min.threshold = 100) {
  # Function to extract 96 mutation type and spectrum channels
  # muts is the data-frame of mutations
  # ls.merge.out is the output from the merge.clusters.differ.by.seq.error() function
  # min.threshold is the minimum number of mutations in cluster for including in heatmap
  
  muts$TYPE <- rep("C>T", nrow(muts))
  muts$TYPE[(muts$ref=="C" & muts$alt == "G") | (muts$ref == "G" & muts$alt == "C")] <- "C>G"
  muts$TYPE[(muts$ref=="C" & muts$alt == "A") | (muts$ref == "G" & muts$alt == "T")] <- "C>A"
  muts$TYPE[(muts$ref=="T" & muts$alt == "A") | (muts$ref == "A" & muts$alt == "T")] <- "T>A"
  muts$TYPE[(muts$ref=="T" & muts$alt == "C") | (muts$ref == "A" & muts$alt == "G")] <- "T>C"
  muts$TYPE[(muts$ref=="T" & muts$alt == "G") | (muts$ref == "A" & muts$alt == "C")] <- "T>G"
  muts$CONTEXT[muts$ref == "G" | muts$ref == "A"] <- 
    sapply(muts$CONTEXT[muts$ref == "G" | muts$ref == "A"], rev.comp)
  muts$spec96 <- paste0(muts$TYPE, ".", substring(muts$CONTEXT, first = 10, last = 12))
  
  muts$cluster <- ls.merge.out[["final.cluster"]]
  cluster.by.spec <- table(muts$cluster, muts$spec96)
  if (plot.heatmap) {
    # temp.clust.spec <- cluster.by.spec[rowSums(cluster.by.spec) >= min.threshold,]
    # row.names(temp.clust.spec) <- paste0("Cl.", row.names(temp.clust.spec))
    idx = which(rowSums(cluster.by.spec) >= min.threshold)
    temp.clust.spec <- t(as.matrix(cluster.by.spec[idx,]))
    if(dim(temp.clust.spec)[1] > dim(temp.clust.spec)[2]){
     temp.clust.spec = t(temp.clust.spec)
    }
    row.names(temp.clust.spec) <- paste0("Cl.", row.names(cluster.by.spec)[idx])

    # to handle the case where there is only 1 row, can't be clustered in heatmap step below
    if(dim(temp.clust.spec)[1] == 1){
     temp.clust.spec = rbind(temp.clust.spec,temp.clust.spec)
    }

    # ensure that there are 96 trinuc contexts, pad with zeros if columns are missing
    if(length(unique(sort(colnames(temp.clust.spec)))) < 96){
     ids.trinuc = read.table(paste(repo_location,'trinuc.contexts.txt',sep='/'))
     ids.trinuc = as.character(ids.trinuc[,1])
     idx = grep(paste(colnames(temp.clust.spec),collapse='|'),ids.trinuc,invert=T)
     for(i in 1:length(idx)){
      temp.clust.spec = cbind(temp.clust.spec,0)
      colnames(temp.clust.spec)[dim(temp.clust.spec)[2]] = ids.trinuc[idx[i]]
     }
     idx = match(ids.trinuc,colnames(temp.clust.spec))
     temp.clust.spec = temp.clust.spec[,idx]
    }

    # condition to handle no clustering scenario due to low similarity in cosine distance
    if(length(as.dist(1 - distance(matrix(c(temp.clust.spec), ncol=96) / rowSums((temp.clust.spec)), method="cosine"))) != 0){

    temp.dendro <- as.dendrogram(hclust(as.dist(1 - distance(matrix(c(temp.clust.spec), ncol=96) / rowSums((temp.clust.spec)), method="cosine"))))

    } else {

     temp.dendro = NULL

    }
    
    heatmap(temp.clust.spec, scale="row",  
            col=brewer.pal(9, "PuBuGn"), Colv=NA, Rowv = temp.dendro, ylab="Cluster",
            xlab = "Mutation type", main = "96 channel mutation spectrum by cluster",
            ColSideColors = rep(brewer.pal(6, "BrBG"), each=16),
            RowSideColors = as.character(cut(rowSums(temp.clust.spec), breaks = 9, 
                                             labels = brewer.pal(9,"Greens"))))
  }
  
  return(muts)
  
}

rev.comp <- function(x) {
  x <- toupper(x)
  x <- rev(unlist(strsplit(x,"")))
  x <- paste(sapply(x, switch,  "A"="T", "T"="A","G"="C","C"="G"),collapse="")
  return(x)
}

#############################################

# Posterior distribution of cluster locations after resolving label switching

post.param.dist_posthoc <- function(gs.out, ls.merge.out, centiles = c(0.025, 0.5, 0.975), ls.type = "ECR", burn.in = 2000, samp.names, plot.heatmap = TRUE, min.threshold = 100) {
  # Function to return centiles of posterior distribution for parameters
  # gs.out is output from Gibbs sampler
  # ls.merge.out is output from the merge.clusters.differ.by.seq.error() function
  # ls.type is the preferred method for label switching correction to be used
  # burn.in is the number of iterations of burn-in supplied to the label.switching() function
  # min.threshold is the minimum number of mutations in a cluster for plotting in heatmap
  
  pi.h.j <- gs.out[["pi.h.j"]]
  num.clusters <- dim(pi.h.j)[3]
  num.samples <- dim(pi.h.j)[2]
  num.iter <- dim(pi.h.j)[1]
  
  pi.mcmc <- ls.merge.out[["pi.mcmc"]]
  ls.final <- ls.merge.out[["final.cluster"]]
  
  result <- array(NA, dim=c(num.clusters, num.samples, length(centiles)),
                  dimnames = list(clusters = paste0("Cl.", 1:num.clusters),
                                  samples = samp.names,
                                  centiles = paste0("Centile.", centiles)))
  
  for (i in 1:num.clusters) {
    for (j in 1:num.samples) {
      result[i,j,] <- quantile(pi.mcmc[,i,j], probs = centiles, na.rm = TRUE)
    }
  }
  
  num.muts <- table(ls.final)
  
  if (plot.heatmap) {
    idx = which(num.muts >= min.threshold)
    if(length(idx) == 1){
     tmp.id = gsub('^','Cl.',names(num.muts)[which(num.muts > min.threshold )])
     temp.r <- t(as.matrix(result[as.double(names(num.muts))[num.muts >=min.threshold],,"Centile.0.5"]))
     rownames(temp.r) = tmp.id
    } else if(length(idx) > 1){
     temp.r <- as.matrix(result[as.double(names(num.muts))[num.muts >=min.threshold],,"Centile.0.5"])
    } else {
     message('No clusters > 100 mutations.')
    }
    # temp.r <- as.matrix(result[as.double(names(num.muts))[num.muts >=min.threshold],,"Centile.0.5"])

    # to handle case of NA clusters in post.param.dist_posthoc
    # but NA clusters should disappear by increasing iterations to 16k
    # idx.na = which(is.na(rowSums(temp.r)))
    # if(length(idx.na) > 0){
    #  temp.r <- temp.r[-idx.na,]
    # }

    # to handle the condition of single row or column matrices that cannot be clustered during heatmap creation
    RowSideColors = as.character(cut( num.muts[num.muts >= min.threshold], breaks = 9, labels = brewer.pal(9,"Greens")))
    if(dim(temp.r)[1] == 1){
     temp.r = rbind(temp.r,temp.r)
    } else if(dim(temp.r)[2] == 1){
     temp.r = cbind(temp.r,temp.r)
     # RowSideColors = rep(RowSideColors,num.samples)
    }

    # condition to handle clusters with N/A mean.vafs there weren't merged: 
    # replace N/A's with zero vafs
    # temp.r[is.na(temp.r)] = 0
    # remove clusters with N/A rows of vafs
    idx = which(rowSums(is.na(temp.r)) > 0)
    if(length(idx) > 0){
     temp.r.bk = temp.r
     temp.r = as.matrix(temp.r[-idx,])
     num.muts = num.muts[num.muts >= min.threshold][-idx]
     if(length(colnames(temp.r)) == 0){
      temp.r = t(temp.r)
      rownames(temp.r) = rownames(temp.r.bk)[-idx]
     }
     if(dim(temp.r)[1] == 1){
      temp.r = rbind(temp.r,temp.r)
     } else if(dim(temp.r)[2] == 1){
      temp.r = cbind(temp.r,temp.r)
     }
    }

    # ensure that the number of colors in RowSideColors equals number of rows (clusters) in temp.r
    # RowSideColors = as.character(cut(seq(1,dim(temp.r)[1]), breaks = 9, labels = brewer.pal(9,"Greens")))

    # condition to handle no clustering scenario due to low similarity in cosine distance
    if(length(as.dist(1 - distance(matrix(c(temp.r), ncol=num.samples) / rowSums((temp.r)), method="cosine"))) != 0){

    temp.dendro <- as.dendrogram(hclust(as.dist(1 - distance(matrix(c(temp.r), ncol=num.samples) / rowSums((temp.r)), method="cosine"))))
    
    } else {

     temp.dendro = NULL

    }

    # heatmap(temp.r, Rowv = temp.dendro, xlab="Sample", ylab="Cluster", main="Posterior median VAFs by cluster and sample", scale="none", col=brewer.pal(9, "Oranges"), RowSideColors = as.character(cut( num.muts[num.muts >= min.threshold], breaks = 9, labels = brewer.pal(9,"Greens"))))
    # heatmap(temp.r, Rowv = temp.dendro, xlab="Sample", ylab="Cluster", main="Posterior median VAFs by cluster and sample", scale="none", col=brewer.pal(9, "Oranges"), RowSideColors = RowSideColors)
    require(gplots)
    heatmap.2(temp.r, Rowv = temp.dendro, ylab="Cluster",
              main="Posterior median VAFs by cluster and sample",
              scale="none", col=brewer.pal(9, "Oranges"),
              RowSideColors = RowSideColors,
              key=T,
              key.xlab = "VAF",
              key.title = "",
              density.info="none",  # turns off density plot inside colour legend
              trace="none",
              margins=c(10,8))
  }
  
  # return(list(centiles = result, num.muts = num.muts, pi.mcmc = pi.mcmc, heat.dat = temp.r, idx.na = idx.na))
  return(list(centiles = result, num.muts = num.muts, pi.mcmc = pi.mcmc, heat.dat = temp.r))
}

#############################################

# Density plots for cluster locations

cluster.density.plot_posthoc <- function(gs.out, ls.merge.out, post.ls.dist, curr.cluster, 
                                 vaf.threshold = 0.05, samp.names) {
  y1 <- gs.out[["y1"]]
  N1 <- gs.out[["N1"]]
  vaf <- y1/N1
  colnames(vaf) <- samp.names
  which.muts <- which(ls.merge.out[["final.cluster"]] == curr.cluster)
  
  pi.mcmc <- post.ls.dist[["pi.mcmc"]][,curr.cluster,]
  colnames(pi.mcmc) <- samp.names
  
  # Find which samples have median VAF greater than threshold; and add one random sample not part of cluster
  which.samples <- which(apply(pi.mcmc, MARGIN=2, FUN = median, na.rm=TRUE) >= vaf.threshold)

  # if (length(which.samples) == 0) {which.samples <- sample(1:ncol(pi.mcmc), size=5)}

  # to handle to case where the number of samples (columns) is < 5
  # if (length(which.samples) == 0 & ncol(pi.mcmc) >= 5){
  #  which.samples <- sample(1:ncol(pi.mcmc), size=5)
  # } else {
  #  which.samples <- sample(1:ncol(pi.mcmc), size=ncol(pi.mcmc))
  # }

  if (length(which.samples) < ncol(pi.mcmc)) {
    # which.samples <- c(which.samples, sample((1:ncol(pi.mcmc))[-which.samples], size = 1))
    which.samples <- c(which.samples, as.numeric(sample(as.character(setdiff(seq(1,dim(pi.mcmc)[2]),which.samples)), size = 1)))
  }
  
  pi.df <- pi.mcmc[, which.samples]
  pi.df <- cbind(pi.df, data.frame(pi.est = factor(rep(20, nrow(pi.df))), raw = factor(rep(NA, nrow(pi.df)))))
  vaf.df <- vaf[which.muts, which.samples]
  vaf.df <- cbind(vaf.df, data.frame(pi.est = factor(rep(NA, nrow(vaf.df))), raw = factor(rep(20, nrow(vaf.df)))))
  combined.df <- rbind(pi.df, vaf.df)
  
  prs.plt <- ggpairs(vaf.df, columns = 1:(ncol(combined.df)-2), 
                     title = paste0("Cluster ", curr.cluster),
                     lower = list(continuous = wrap("points", alpha = 0.3, size=0.25, col="plum4"))) +
                     #upper = list(continuous = wrap("density"))) +
    theme(panel.grid.minor = element_blank()) 
  pm2 <- prs.plt
  for(i in 2:prs.plt$nrow) {
    for(j in 1:(i-1)) {
      pm2[i,j] <- prs.plt[i,j] +
        scale_x_continuous(limits = c(0, 0.75)) +
        scale_y_continuous(limits = c(0, 0.75))
    }
  }
  
  return(pm2)
}



