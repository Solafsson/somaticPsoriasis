rev.comp <- function(x)
{
  x <- toupper(x)
  x <- rev(unlist(strsplit(x, "")))
  x <- paste(sapply(x, switch, A = "T", T = "A", G = "C", 
                    C = "G"), collapse = "")
  return(x)
}

context.extract_posthoc <- function(muts, ls.merge.out, plot.heatmap = T, min.threshold = 15) {
  # Function to extract 96 mutation type and spectrum channels
  # muts is the data-frame of mutations
  # ls.merge.out is the output from the merge.clusters.differ.by.seq.error() function
  # min.threshold is the minimum number of mutations in cluster for including in heatmap
  message("Inside context.extract_posthoc")
  print("Type defined")
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
  print("Here I am!")
  if (plot.heatmap) {
    # temp.clust.spec <- cluster.by.spec[rowSums(cluster.by.spec) >= min.threshold,]
    # row.names(temp.clust.spec) <- paste0("Cl.", row.names(temp.clust.spec))
    idx = which(rowSums(cluster.by.spec) >= min.threshold)
    temp.clust.spec <- t(as.matrix(cluster.by.spec[idx,]))
    if(dim(temp.clust.spec)[1] > dim(temp.clust.spec)[2]){
      temp.clust.spec = t(temp.clust.spec)
    }
    #row.names(temp.clust.spec) <- paste0("Cl.", row.names(cluster.by.spec)[idx])
    rownames(temp.clust.spec) <- paste("Cl.", rownames(temp.clust.spec), sep="")
    
    # to handle the case where there is only 1 row, can't be clustered in heatmap step below
    if(dim(temp.clust.spec)[1] == 1){
      temp.clust.spec = rbind(temp.clust.spec,temp.clust.spec)
    }
    
    # ensure that there are 96 trinuc contexts, pad with zeros if columns are missing
    if(length(unique(sort(colnames(temp.clust.spec)))) < 96){
      ids.trinuc = read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/06_phylogenics/hdp_clustering/snv-clustering-using-the-dirichlet-process/trinuc.contexts.txt")
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


### get mutational spectrum plots
extract.mut.spec <- function(wd) {
  files = list.files(path=wd,pattern='Rsession_ls2.dat',full.names=T,recursive=T)
  for(n in 1:length(files)) {
    load(files[n])
    wd=outdir=getwd()
    mut.spec.by.clust <- context.extract_posthoc(muts=mut.spec, ls.merge.out=lymph.merge.ls,
                                                 plot.heatmap = T, min.threshold = 15)
    #save(mut.spec.by.clust,file=paste(gsub('\\/Rsession_ls2.dat','',files[n]),'mut.spec.by.clust.RData',sep='/'))
  }
  return(mut.spec.by.clust)
}



### assign mutations to clusters
trunc <- function(x, ..., prec = 0) base::trunc(x * 10^prec, ...) / 10^prec

# append new_cluster_id to each mutation per cluster
cluster.append <- function(snvdir, ndpdir, mut.spec.by.clust,id.patient) {
  x.cluster_and_samples = fread(paste0(ndpdir,'/cluster_and_samples.csv'))
  rownames(mut.spec.by.clust) = paste(mut.spec.by.clust$chrom,mut.spec.by.clust$pos,mut.spec.by.clust$ref,mut.spec.by.clust$alt,sep='_')
  mut.spec.by.clust$cluster = gsub('^','Cl.',mut.spec.by.clust$cluster)
  snv_tbl = NULL
  idx = grep('P',colnames(x.cluster_and_samples))
  ids = colnames(x.cluster_and_samples)[idx]
  
  file.snv = paste(snvdir,paste0(id.patient,'_bb_pass_snvs_all.csv'),sep='/')
  if(!file.exists(file.snv)){next}
  tmp = read.table(file.snv,header=T,sep=',')
  tmp$cluster_id = NA
  for(l in 1:dim(mut.spec.by.clust)[1]){
    idx = which(tmp$mut_id %in% rownames(mut.spec.by.clust)[l])
    tmp$cluster_id[idx] = mut.spec.by.clust$cluster[l]
  }
  snv_tbl = rbind(snv_tbl,tmp)

  return(x.cluster_and_samples)
}

# save num.muts
extract.mutation.counts <- function(x.cluster_and_samples) {
  num.muts = x.cluster_and_samples$no.of.mutations.assigned
  names(num.muts) = x.cluster_and_samples$cluster_id
  num_df <- as.data.frame(num.muts)
  num_df <- rownames_to_column(num_df, var = "cluster_id")
  return(num_df)
} 

# extract 95% CI VAF data
extract.95.CI <- function(x.cluster_and_samples, ndpdir) {
  x.clusters.unique.fullrun.centiles = NULL
  colnames(x.cluster_and_samples) = gsub('^PD.*._lo','lo',colnames(x.cluster_and_samples))
  dirichlet_obj = readRDS(file.path(paste(ndpdir,'object_post.cluster.pos.dat',sep='/')))
  for(j in 1:dim(x.cluster_and_samples)[1]){
    tmp = as.data.frame(dirichlet_obj$centiles[as.numeric(as.character(gsub('Cl.','',x.cluster_and_samples$cluster_id)))[j],,])
    tmp$id = gsub('^PD.*._lo','lo',rownames(tmp))
    tmp$no.of.mutations.assigned = x.cluster_and_samples$no.of.mutations.assigned[j]
    tmp$cluster_id = x.cluster_and_samples$cluster_id[j]
    x.clusters.unique.fullrun.centiles = rbind(x.clusters.unique.fullrun.centiles,tmp)
  }
  rownames(x.clusters.unique.fullrun.centiles) = NULL
  return(x.clusters.unique.fullrun.centiles)
}

# prep median vaf table 
create.med.vaf.tbl <- function(x.cluster_and_samples) {
  x.cluster_and_samples = as.matrix(x.cluster_and_samples)
  rownames(x.cluster_and_samples) = x.cluster_and_samples[,'cluster_id']
  idx = grep('P',colnames(x.cluster_and_samples),invert=T)
  x.cluster_and_samples = x.cluster_and_samples[,-idx]
  mode(x.cluster_and_samples) = 'numeric'
  x.cluster_and_samples = trunc(x.cluster_and_samples, prec = 2)
  return(x.cluster_and_samples)
}

# pigeonhole plots
plot_cluster_vaf <- function(x.clusters.unique.real.merged, x.clusters.unique.real.centiles.final, target_clust, outdir, id.patient){
  # Generate all cluster pairs
  clust_pairs = get_cluster_pairs(x.clusters.unique.real.merged)
  
  grobs = list()
  centiles_tbl.out = list(data.frame())
  for(i in 1:dim(clust_pairs)[1]){
    this_pair = clust_pairs[i,]
    cl1_id = this_pair$cl1
    cl2_id = this_pair$cl2
    
    message(sprintf('Plotting cluster pair %s and %s', cl1_id, cl2_id))
    
    target_samp = colnames(x.clusters.unique.real.merged)
    
    cluster_centiles = list()
    for(k in c(1:length(target_samp))) {
      idx1 = which(x.clusters.unique.real.centiles.final$id == ifelse(nchar(target_samp[k]) == 15, str_sub(target_samp[k], -6, -1), target_samp[k])   & x.clusters.unique.real.centiles.final$cluster_id == cl1_id)
      idx2 = which(x.clusters.unique.real.centiles.final$id == ifelse(nchar(target_samp[k]) == 15, str_sub(target_samp[k], -6, -1), target_samp[k]) & x.clusters.unique.real.centiles.final$cluster_id == cl2_id)
      if(length(idx1) == 0 || length(idx2) == 0){
        cl1 = rep(NA,3)
        cl2 = rep(NA,3)
        names(cl1) = c('Centile.0.025','Centile.0.5','Centile.0.975')
        names(cl2) = c('Centile.0.025','Centile.0.5','Centile.0.975')
      } else {
        cl1 = x.clusters.unique.real.centiles.final[idx1,c('Centile.0.025','Centile.0.5','Centile.0.975')]
        cl2 = x.clusters.unique.real.centiles.final[idx2,c('Centile.0.025','Centile.0.5','Centile.0.975')]
      }
      
      cluster_centiles[[k]] = data.frame(
        cl1 = this_pair[1],
        cl2 = this_pair[2],
        sample = target_samp[k],
        cl1_Centile.0.025 = as.numeric(cl1['Centile.0.025']),
        cl1_Centile.0.5 = as.numeric(cl1['Centile.0.5']),
        cl1_Centile.0.975 = as.numeric(cl1['Centile.0.975']),
        cl2_Centile.0.025 = as.numeric(cl2['Centile.0.025']),
        cl2_Centile.0.5 = as.numeric(cl2['Centile.0.5']),
        cl2_Centile.0.975 = as.numeric(cl2['Centile.0.975'])
      )
    }
    
    centiles_tbl = tbl_df(as.data.frame((do.call('rbind', cluster_centiles))))
    centiles_tbl.out[[i]] = as.data.frame(centiles_tbl)
    names(centiles_tbl.out)[i] = paste(cl1_id,cl2_id,sep='_')
    
    grobs[[i]] = centiles_tbl %>% 
      ggplot(aes(x=2*cl1_Centile.0.5, y=2*cl2_Centile.0.5)) +
      geom_errorbar(aes(ymin=2*cl2_Centile.0.025, ymax=2*cl2_Centile.0.975), color='red', width=0) +
      geom_errorbarh(aes(xmin=2*cl1_Centile.0.025, xmax=2*cl1_Centile.0.975), color='red', height=0) +
      theme_light() +
      geom_hline(yintercept = 0.5) +
      geom_vline(xintercept = 0.5) +
      geom_abline(intercept = 0, slope=1, size=0.2) +
      geom_abline(intercept = c(0.7,1), slope=-1, size=0.2, color='grey') +
      # geom_text(aes(label = sample), nudge_x = 0.07, nudge_y = 0.025, size=1) +
      geom_text(aes(label = sample), position=position_jitter(width=0.025,height=0.025), size=1) +
      scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
      scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
      labs(x=sprintf('Cell fraction, cluster %s.', cl1_id), y=sprintf('Cell fraction, cluster %s.', cl2_id))
     ggsave(file.path(outdir, sprintf('pigeonhole_clusters_%s_%s.pdf', cl1_id, cl2_id)), width=12, height=12, units='cm')
  }
  # Plot all cluster pairs in a single figure
  pdf(file.path(outdir, sprintf('pigeonhole_clusterpairs_%s.pdf', id.patient)), w=20, h=20)
  do.call('grid.arrange', grobs)
  dev.off()
  return(centiles_tbl.out)
}

# determine strength of evidence for inferred clonal relationships
evidence <- function(cl1.vaf, cl2.vaf, flag.nested){
  if(all((cl1.vaf + cl2.vaf) > 1) & flag.nested == T){
    strength = 'strong'
  } else if(all((cl1.vaf + cl2.vaf) >= 0.7) & flag.nested == T){
    strength = 'medium'
  } else if(all((cl1.vaf + cl2.vaf) >= 0.7) & flag.nested == F){
    strength = 'weak'
  } else if(all((cl1.vaf + cl2.vaf) < 0.7)){
    strength = 'very weak'
  } else {
    strength = 'very weak'
  }
  return(strength)
}

# decide which branches to keep
keep <- function(x.branches){
  g = graph_from_data_frame(x.branches,directed=T)
  x.keep = NULL
  idx.letter = 1
  for(j in 1:dim(x.branches)[1]){
    g.paths = all_simple_paths(g,x.branches[j,1],x.branches[j,2])
    idx.longest.path = which.max(rapply(g.paths, function(x) length(x)))
    for(i in 1:(length(g.paths[[idx.longest.path]])-1)){
      x.keep = rbind(x.keep,c(names(g.paths[[idx.longest.path]][i]),names(g.paths[[idx.longest.path]][i+1])))
    }
    idx.letter = idx.letter + 1
  }
  colnames(x.keep) = c('From','To')
  x.keep = as.data.frame(x.keep)
  return(x.keep)
}

# create diagnostic plots and gather data
create.diag.plots <- function(x.cluster_and_samples, x.clusters.centiles, outdir, id.patient, num.muts, min_num_mut) {
  target_clust = which(num.muts > min_num_mut)
  out = plot_cluster_vaf(x.cluster_and_samples, x.clusters.centiles, target_clust, outdir, id.patient)
  centiles_tbl_all = na.omit(do.call('rbind', out))
  return(centiles_tbl_all )
}

# determine nesting branches
infer.tree.branches <- function(clust_pairs, centiles_tbl_all, min_contrib = 0.10) {
  x.branches = NULL
  min_contrib = min_contrib
  idx.letter = 1
  for(i in 1:dim(clust_pairs)[1]){
    idx = which(centiles_tbl_all$cl1 == as.character(clust_pairs[i,'cl1']) & centiles_tbl_all$cl2 == as.character(clust_pairs[i,'cl2']) & (centiles_tbl_all$cl1_Centile.0.5 > min_contrib | centiles_tbl_all$cl2_Centile.0.5 > min_contrib))
    idx.cl1 = grep('cl1_Centile',colnames(centiles_tbl_all))
    idx.cl2 = grep('cl2_Centile',colnames(centiles_tbl_all))
    flag.nested = F
    if(all(sign(centiles_tbl_all[idx,idx.cl1] - centiles_tbl_all[idx,idx.cl2]) == 1)){
      # nesting condition 1 detected
      flag.nested = T
      # subset to only include samples present in both clusters
      idx_pres = which(centiles_tbl_all$cl1 == as.character(clust_pairs[i,'cl1']) & centiles_tbl_all$cl2 == as.character(clust_pairs[i,'cl2']) & (centiles_tbl_all$cl1_Centile.0.5 > min_contrib & centiles_tbl_all$cl2_Centile.0.5 > min_contrib))
      if (length(idx_pres) > 0){
        x.branches = rbind(x.branches,c(unique(sort(centiles_tbl_all$cl1[idx_pres])),unique(sort(centiles_tbl_all$cl2[idx_pres])),evidence(centiles_tbl_all[idx_pres,idx.cl1][2], centiles_tbl_all[idx_pres,idx.cl2][2], flag.nested)))
      }
    } else if(!all(sign(centiles_tbl_all[idx,idx.cl1] - centiles_tbl_all[idx,idx.cl2]) == 1) & all(centiles_tbl_all[idx,idx.cl1][2] - centiles_tbl_all[idx,idx.cl2][2] >= -0.03) ){
      flag.nested = F
      # subset to only include samples present in both clusters
      idx_pres = which(centiles_tbl_all$cl1 == as.character(clust_pairs[i,'cl1']) & centiles_tbl_all$cl2 == as.character(clust_pairs[i,'cl2']) & (centiles_tbl_all$cl1_Centile.0.5 > min_contrib & centiles_tbl_all$cl2_Centile.0.5 > min_contrib))
      if (length(idx_pres) > 0){
        x.branches = rbind(x.branches,c(unique(sort(centiles_tbl_all$cl1[idx_pres])),unique(sort(centiles_tbl_all$cl2[idx_pres])),evidence(centiles_tbl_all[idx_pres,idx.cl1], centiles_tbl_all[idx_pres,idx.cl2], flag.nested)))
      }
    } else if(all(sign(centiles_tbl_all[idx,idx.cl1] - centiles_tbl_all[idx,idx.cl2]) == -1)){
      # nesting condition 2 detected
      flag.nested = T
      # subset to only include samples present in both clusters
      idx_pres = which(centiles_tbl_all$cl1 == as.character(clust_pairs[i,'cl1']) & centiles_tbl_all$cl2 == as.character(clust_pairs[i,'cl2']) & (centiles_tbl_all$cl1_Centile.0.5 > min_contrib & centiles_tbl_all$cl2_Centile.0.5 > min_contrib))
      if (length(idx_pres) > 0){
        x.branches = rbind(x.branches,c(unique(sort(centiles_tbl_all$cl2[idx_pres])),unique(sort(centiles_tbl_all$cl1[idx_pres])),evidence(centiles_tbl_all[idx_pres,idx.cl1][2], centiles_tbl_all[idx_pres,idx.cl2][2], flag.nested)))
      }
    } else if(!all(sign(centiles_tbl_all[idx,idx.cl1] - centiles_tbl_all[idx,idx.cl2]) == -1) & all(centiles_tbl_all[idx,idx.cl1][2] - centiles_tbl_all[idx,idx.cl2][2] <= 0.03)){
      flag.nested = F
      # subset to only include samples present in both clusters
      idx_pres = which(centiles_tbl_all$cl1 == as.character(clust_pairs[i,'cl1']) & centiles_tbl_all$cl2 == as.character(clust_pairs[i,'cl2']) & (centiles_tbl_all$cl1_Centile.0.5 > min_contrib & centiles_tbl_all$cl2_Centile.0.5 > min_contrib))
      if (length(idx_pres) > 0){
        x.branches = rbind(x.branches,c(unique(sort(centiles_tbl_all$cl2[idx_pres])),unique(sort(centiles_tbl_all$cl1[idx_pres])),evidence(centiles_tbl_all[idx_pres,idx.cl1], centiles_tbl_all[idx_pres,idx.cl2], flag.nested)))
      }
    } else {
      # no nesting condition detected
      x.branches = rbind(x.branches,c('P',unique(sort(centiles_tbl_all$cl1[idx])),'strong'))
      x.branches = rbind(x.branches,c('P',unique(sort(centiles_tbl_all$cl2[idx])),'strong'))
    }
  }
  colnames(x.branches) = c('From','To','Evidence')
  x.branches = unique(x.branches)
  
  # append independent clusters
  idx = grep('P',unique(x.branches[,'From']),invert=T)
  for(i in 1:length(idx)){
    x.branches = rbind(x.branches,c('P',as.character(unique(x.branches[,'From'])[idx[i]]),'strong'))
  }
  
  # append indep clusters not captured in clust_pairs
  idx = grep(paste(unique(x.branches[,'From']),collapse='|'),rownames(x.cluster_and_samples.final),invert=T)
  for(i in 1:length(idx)){
    x.branches = rbind(x.branches,c('P',rownames(x.cluster_and_samples.final)[idx[i]],'strong'))
  }
  x.branches = unique(x.branches)
  
  return(x.branches)
}
 
# break ties if multiple potential paths, going with highest med vaf
break.ties <- function(x.branches.keep, centiles_tbl_all, min_vaf_threshold = 0.10) {
  branch_count <- table(unlist(x.branches.keep[,2]))
  multi_branches <- names(branch_count[branch_count>1])
  message(multi_branches, " have multiple possible paths")
  if (length(multi_branches>0)) {
    for (i in 1:length(multi_branches)) {
      problem_parents <- x.branches.keep[x.branches.keep$To == multi_branches[i],]$From
      med_vaf <- list()
      for (j in 1:length(problem_parents)) {
        clust_vafs <- centiles_tbl_all[centiles_tbl_all$cl2==problem_parents[j] & centiles_tbl_all$cl1==multi_branches[i],8]
        if (length(clust_vafs) >0) {
          idx_pres <- which(clust_vafs >min_vaf_threshold )
          med_vaf[j] <- median(clust_vafs[idx_pres])
        } else {
          clust_vafs <- centiles_tbl_all[centiles_tbl_all$cl1==problem_parents[j] & centiles_tbl_all$cl2==multi_branches[i],5]
          idx_pres <- which(clust_vafs >min_vaf_threshold )
          med_vaf[j] <- median(clust_vafs[idx_pres])
        }
      }
      winning_idx <- which.max(med_vaf)
      winning_parent <- problem_parents[winning_idx]
      losing_parents <- problem_parents[!problem_parents %in% winning_parent]
      x.branches.keep <- x.branches.keep %>%
        dplyr::filter(!(From %in% losing_parents & To == multi_branches[i]))
    } 
  }
  return(x.branches.keep)
}
