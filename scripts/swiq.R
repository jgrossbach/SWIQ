swiq <- function(y,
                 x,
                 num.trees=10000,
                 num.threads=1,
                 min.bucket=15,
                 imp.blocks=1,
                 max.trees=Inf,
                 min.trees=0,
                 pop.str=NULL,
                 pick.threshold=1){
  #load requirements
  require(ranger)
  require(abind)
  require(parallel)
  
  # check all input
  # if(is.null(num.threads)){
  #   num.threads <- imp.blocks
  # }
  
  # preprocess data
  if(!is.data.frame(x)){
    x <- data.frame(x)
  }
  y <- scale(y,center = T,scale = T)[,1]
  if(any(is.na(y))){
    na_bol <- is.na(y)
    y <- y[!na_bol]
    x <- x[!na_bol,]
  }
  
  if(length(y)<(2*min.bucket)){
    return(rep(1,ncol(x)))
  }
  # in blocks: determine the effect on accuracy of each marker for each sample 
  # and aggregate
  trees_by_block <- num.trees
  
  imp_mat_list <- lapply(1:imp.blocks,FUN=function(i){
    #i <<- i
    get_importances(x=x,
                    y=y,
                    num.trees.block=trees_by_block,
                    num.threads=num.threads,
                    min.bucket=min.bucket,
                    max.trees=max.trees,
                    min.trees=min.trees,
                    pop.str=pop.str,
                    pick.threshold=pick.threshold)
  })
  
  # aggregate blocks
  imps <- Reduce("+",imp_mat_list)
  
  # test each marker
  n_pos <- rowSums(imps>0,na.rm=T)
  n <- ncol(imps)
  # only test if the number of samples with improved predictions is bigger
  # than expected, not smaller
  p <- pbinom(n_pos-1,size = n,prob = 0.5,lower.tail = F)
  # p <- sapply(1:nrow(imps),FUN=function(m){
  #   m <<- m
  #   global_alleles <- as.numeric(na.omit(unique(x[,m])))
  #   pvs <- sapply(global_alleles,FUN=function(a){
  #     binom.test(sum(imps[m,which(x[,m]==a)]>0),
  #                n = length(which(x[,m]==a)),
  #                alternative = "greater")$p.value
  #   })
  #   p_m <- 1-prod(1-pvs)
  #   # p_m <- fisher_method(pvs)
  #   return(p_m)
  # })
  return(p)
}

get_importances <- function(x,
                            y,
                            num.trees.block,
                            num.threads,
                            min.bucket,
                            max.trees,
                            min.trees,
                            pop.str,
                            pick.threshold){
  #impute NAs, including pop structure predictors
  x <- apply(x,2,FUN=function(m){
    if(any(is.na(m))){
      m[is.na(m)] <- sample(m[!is.na(m)],replace = T,size = sum(is.na(m)))
    }
    return(m)
  })
  x <- as.data.frame(x)
  nx <- ncol(x) # number of markers (not pop structure)
  if(!is.null(pop.str)){
    x <- cbind(x,pop.str)
    pop.str.ind <- (nx+1):ncol(x)
    pop.str.names <- colnames(x)[pop.str.ind]
  }else{
    pop.str.ind <- NULL
    pop.str.names <- NULL
  }
  
  
  
  # make the forest
  rf <- ranger(y = y,
               x = x,
               importance="none",
               num.trees = num.trees.block,
               num.threads = num.threads,
               min.bucket = min.bucket,
               keep.inbag = T,
               always.split.variables = pop.str.names)
  
  # compute in which trees a sample is oob
  # this has to be computed before deforestation
  oob <- do.call("cbind",rf$inbag.counts)==0
  
  # which markers were used where?
  m_use_orig <- mclapply(1:rf$num.trees,mc.cores = num.threads,FUN=function(nt){
    (1:ncol(x))%in%(treeInfo(rf,nt)$splitvarID+1)
  })
  m_use_orig <- do.call("cbind",m_use_orig)
  
  # Do not consider population structure predictors
  # for which trees are interesting
  m_use_orig[pop.str.ind,] <- F
  
  # how often was each marker used?
  n_picked <- rowSums(m_use_orig[1:nx,])
  med_picked <- median(n_picked)
  
  # which markers were used more often than the median?
  p_picked <- pbinom(q = n_picked-1,
                     size = rf$num.trees,
                     prob = med_picked/rf$num.trees,
                     lower.tail = F)
  picked <- p_picked<=pick.threshold
  
  # ignore the markers that were picked not often enough
  m_use_orig[1:nx,][!picked,] <- F
  
  
  
  # reduce the forest to the usefull trees
  m_use_cum <- t(apply(m_use_orig,1,cumsum))
  m_use_max <- m_use_cum[,rf$num.trees]
  if(any(m_use_max<min.trees)){
    m_use_cum[m_use_max<min.trees,] <- 0
    m_use_orig[m_use_max<min.trees,] <- F
  }
  usefull <- c(T,sapply(2:rf$num.trees,FUN=function(nt){
    any((m_use_cum[,nt]>m_use_cum[,nt-1])&m_use_cum[,nt]<(max.trees+1))
  }))
  rf <- suppressWarnings(deforest(rf,which(!usefull)))
  m_use <- m_use_orig[,usefull]
  oob <- oob[,usefull]
  
  # identify the first x trees each marker was used
  too_many_trees <- rowSums(m_use)>max.trees
  if(any(too_many_trees)){
    for(m in which(too_many_trees)){
      last_tree <- which(m_use[m,])[max.trees]
      m_use[m,(last_tree+1):rf$num.trees] <- F
    }
  }
  
  # compute in which node each sample would end up in
  preds_real <- predict(object = rf,
                        data = x,
                        predict.all = T,
                        type="terminalNodes",
                        num.threads=num.threads)$predictions
  max_node <- max(preds_real)
  empty_out <- rep(0,max_node)
  
  
  # extract the prediction of each terminal node
  node_pred <- mclapply(1:rf$num.trees,mc.cores = num.threads,FUN=function(nt){
    out <- rep(NA,max_node)
    tI <- treeInfo(rf,nt)
    tI <- tI[tI$terminal,]
    out[tI$nodeID] <- tI$prediction
    return(out)
  })
  node_pred <- do.call("cbind",node_pred)
  
  # compute the raw error for each node using all OOB samples that end up in 
  # this node without changing any alleles
  node_error_means <- mclapply(1:rf$num.trees,
                               mc.cores = num.threads,
                               FUN=function(nt){
                                 out <- rep(NA,max_node)
                                 tI <- treeInfo(rf,nt)
                                 tI <- tI[tI$terminal,]
                                 for(i in 1:nrow(tI)){
                                   n <- tI[i,"nodeID"]
                                   # TODO break into two lines
                                   res <- (y-tI[i,"prediction"])[(preds_real[,nt]==n)&oob[,nt]]
                                   err_adj <- (y-mean(y))[(preds_real[,nt]==n)&oob[,nt]]
                                   out[n] <- mean(res^2)
                                 }
                                 return(out)
                               })
  node_error_means <- do.call("cbind",node_error_means)
  
  # TODO check if necessary, do i ever use the IB node sizes?
  node_size <- mclapply(1:rf$num.trees,mc.cores = num.threads,FUN=function(nt){
    out <- rep(NA,max_node)
    tI <- treeInfo(rf,nt)
    tI <- tI[tI$terminal,]
    for(i in 1:nrow(tI)){
      n <- tI[i,"nodeID"]
      ct <- sum(preds_real[,nt]==n&!oob[,nt])
      out[n] <- ct
    }
    return(out)
  })
  node_size <- do.call("cbind",node_size)
  
  # compute how many OOB samples end up in each node, without altering alleles
  node_size_oob <- mclapply(1:rf$num.trees,mc.cores = num.threads,FUN=function(nt){
    out <- rep(NA,max_node)
    tI <- treeInfo(rf,nt)
    tI <- tI[tI$terminal,]
    for(i in 1:nrow(tI)){
      n <- tI[i,"nodeID"]
      ct <- sum(preds_real[,nt]==n&oob[,nt])
      out[n] <- ct
    }
    return(out)
  })
  node_size_oob <- do.call("cbind",node_size_oob)
  
  # compute the frequencies by tree and allele for OOB
  global_alleles <- sort(unique(unlist(x[,setdiff(1:ncol(x),pop.str.ind)])))
  amat <- lapply(global_alleles,FUN=function(a){
    x[,setdiff(1:ncol(x),pop.str.ind)]==a
  })
  amat <- abind(amat,along=3)
  
  # TODO can i get rid of the AFs completely?
  tx <- t(x)
  af_by_tree <- mclapply(1:rf$num.trees,mc.cores = num.threads,FUN=function(nt){
    nt <<- nt
    x_t <- t(tx[,oob[,nt]])
    allele_freqs <- sapply(1:length(global_alleles),FUN=function(i){
      colMeans(amat[,,i][oob[,nt],])
    })
    return(allele_freqs)
  })
  af_by_tree <- abind(af_by_tree,along = 3)
  
  
  # compute how important each marker is for the prediction of each sample
  imps <- mclapply(1:ncol(x),mc.cores = num.threads,FUN = function(m){
    #m <<- m
    use_trees <- m_use[m,]
    
    # markers that were not picked often enough were set to F and can be 
    # filtered out here
    if(sum(use_trees)==0){
      return(rep(0,nrow(x)))
    }else{
      
      # keep only trees where the marker was used
      rf_sub <- suppressWarnings(deforest(rf,!use_trees))
      
      # extract the AFs for this marker across the trees
      afreq_m <- af_by_tree[m,,]
      
      # make a copy of the genotype
      x_m <- x
      
      # compute in which node a sample would end up in given each 
      # possible allele
      node_array <- lapply(global_alleles,FUN=function(a){
        x_m[,m] <- a
        nodes_a <- predict(object = rf_sub,
                           data = x_m,
                           predict.all = T,
                           type="terminalNodes",
                           num.threads=1)$predictions
        return(nodes_a)
      })
      node_array <- abind(node_array,along=3)
      
      # for each sample extract the OOB node sizes of the node where the sample
      # ends up in for each allele
      # for the true node of the sample, the size is reduced by one to count
      # only other samples
      node_pops <- lapply(1:length(global_alleles),FUN=function(i){
        pos_mat <- cbind(as.vector(node_array[,,i]),
                         rep(which(use_trees),each=length(y)))
        pops <- matrix(node_size_oob[pos_mat],ncol = ncol(node_array))
        true_allele <- node_array[,,i]==preds_real[,use_trees]
        pops[true_allele] <- pops[true_allele]#-1
        return(pops)
      })
      allele_sum <- Reduce("+",node_pops)
      node_pops <- lapply(node_pops,FUN=function(x){
        x/allele_sum
      })
      node_pops <- abind(node_pops,along=3)
      
      preds_real_sub <- preds_real[,use_trees]
      a_preds <- lapply(1:length(global_alleles),FUN=function(i){
        a <- global_alleles[i]
        if(sum(afreq_m[i,])==0){
          return(matrix(0,nrow=nrow(x),ncol=rf_sub$num.trees))
        }
        x_m[,m] <- a
        nodes_a <- predict(object = rf_sub,
                           data = x_m,
                           predict.all = T,
                           type="terminalNodes",
                           num.threads=1)$predictions
        pos_mat <- cbind(as.vector(nodes_a),
                         rep(which(use_trees),each=length(y)))
        y_a <- matrix(node_pred[pos_mat],
                      nrow=length(y))
        error_estimate_a <- matrix(node_error_means[pos_mat],
                                   nrow=length(y))
        #compute an adjustment for the error by computing oob errors without
        # the sample of interest
        size_mat <- matrix(node_size_oob[pos_mat],
                           nrow=length(y))
        adj_error_a <- (error_estimate_a*size_mat-(y_a-y)^2)/(size_mat-1)
        # adj_error_a <- error_estimate_a-adj_error_a
        oob_a <- matrix(oob[,use_trees],
                        nrow=length(y))
        oob_a <- oob_a&(nodes_a==preds_real_sub)
        # error_estimate_a <- (error_estimate_a*!oob_a)+(adj_error_a*oob_a)
        # error_estimate_a[is.infinite(error_estimate_a)] <- NA
        #error_estimate_a[oob_a] <- adj_error_a[oob_a]
        error_estimate_a[is.infinite(error_estimate_a)] <- NA
        error_a <- (y_a-y)^2
        # oob_m <- oob[,use_trees,drop=F]
        # empty_in_node <- matrix(F,nrow = nrow(x),ncol = max_node)
        # for(nt in 1:rf_sub$num.trees){
        #   
        #   hit_mat <- cbind(1:nrow(x),nodes_a[,nt])
        #   hit_mat <- hit_mat[oob_m[,nt],]
        #   
        #   error_vec <- error_a[,nt]
        #   node_error_means_nt <- node_error_means[,use_trees][,nt]
        #   # for(j in 1:nrow(hit_mat)){
        #   #   node_error_means[hit_mat[j,2]] <- node_error_means[hit_mat[j,2]]+error_vec[hit_mat[j,1]]
        #   # }
        #   # node_counts <- tabulate(hit_mat[,2],max_node)
        #   # node_error_means <- node_error_means/node_counts
        #   
        #   error_a[,nt] <- error_a[,nt]-node_error_means_nt[nodes_a[,nt]]
        # }
        
        # return(error_a-error_estimate_a)
        return(error_a)
      })
      
      a_preds <- abind(a_preds,along = 3)
      
      real_error <- matrix(NA,nrow = nrow(x),ncol = rf_sub$num.trees)
      blended_error <- matrix(0,nrow = nrow(x),ncol = rf_sub$num.trees)
      for(i in 1:length(global_alleles)){
        use_samples <- x[,m]==global_alleles[i]
        real_error[use_samples,] <- a_preds[use_samples,,i]
        
        blended_error <- blended_error+a_preds[,,i]*node_pops[,,i]
      }
      diff_error <- blended_error-real_error
      diff_error[!oob[,use_trees,drop=F]] <- 0
      #some errors are NA because there were no oob samples in the node
      diff_error[is.na(diff_error)] <- 0
      same <- apply(a_preds,
                    seq_len(length(dim(a_preds)) - 1),
                    function(z) length(unique(z)) == 1)
      diff_error[same] <- 0 #remove numerical issues
      
      # #weight for which samples are part of the training data to punish outliers
      # sample_use <- rowMeans(oob[,use_trees])/mean(oob)
      # sample_use_mat <- oob*sample_use
      # sample_use_mat[sample_use_mat==0] <- NA
      # sample_corr_fac <- apply(sample_use_mat,2,prod,na.rm=T)[use_trees]
      # diff_error <- t(t(diff_error)/sample_corr_fac)
      
      return(rowMeans(diff_error)*mean(m_use_orig[m,]))
    }
  })
  imps <- do.call("rbind",imps)
  imps <- imps[1:nx,]
  # return the result matrices
  return(imps)
}
