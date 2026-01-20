swiq <- function(y,
                 x,
                 num.trees=10000,
                 num.threads=1,
                 min.bucket=15,
                 imp.blocks=1,
                 max.trees=Inf,
                 min.trees=0){
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
  
  # in blocks: determine the effect on accuracy of each marker for each sample and aggregate
  # can be one block. the total number of trees is distributed
  trees_by_block <- num.trees
  
  imp_mat_list <- lapply(1:imp.blocks,FUN=function(i){
    #i <<- i
    get_importances(x=x,
                    y=y,
                    num.trees.block=trees_by_block,
                    num.threads=num.threads,
                    min.bucket=min.bucket,
                    max.trees=max.trees,
                    min.trees=min.trees)
  })
  
  # aggregate blocks
  imps <- Reduce("+",imp_mat_list)
  
  # test each marker
  n_pos <- rowSums(imps>0,na.rm=T)
  n <- ncol(imps)
  p <- pbinom(n_pos-1,size = n,prob = 0.5,lower.tail = F)
  return(p)
}

get_importances <- function(x,
                            y,
                            num.trees.block,
                            num.threads,
                            min.bucket,
                            max.trees,
                            min.trees){
  #impute NAs
  x <- apply(x,2,FUN=function(m){
    if(any(is.na(m))){
      m[is.na(m)] <- sample(m[!is.na(m)],replace = T,size = sum(is.na(m)))
    }
    return(m)
  })
  x <- as.data.frame(x)
  
  # make the forest
  rf <- ranger(y = y,
               x = x,
               importance="none",
               num.trees = num.trees.block,
               num.threads = num.threads,
               min.bucket = min.bucket,
               keep.inbag = T)
  
  #compute in which trees a sample is oob
  # this has to be computed before deforestation
  oob <- do.call("cbind",rf$inbag.counts)==0
  
  # which markers were used where?
  m_use_orig <- mclapply(1:rf$num.trees,mc.cores = num.threads,FUN=function(nt){
    (1:ncol(x))%in%(treeInfo(rf,nt)$splitvarID+1)
  })
  m_use_orig <- do.call("cbind",m_use_orig)
  
  
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
  
  
  
  # compute the frequencies by tree and allele for OOB
  global_alleles <- sort(unique(unlist(x)))
  amat <- lapply(global_alleles,FUN=function(a){
    x==a
  })
  amat <- abind(amat,along=3)
  
  tx <- t(x)
  af_by_tree <- mclapply(1:rf$num.trees,mc.cores = num.threads,FUN=function(nt){
    #nt <<- nt
    x_t <- t(tx[,oob[,nt]])
    allele_freqs <- sapply(1:length(global_alleles),FUN=function(i){
      colMeans(amat[,,i][oob[,nt],])
    })
    return(allele_freqs)
  })
  af_by_tree <- abind(af_by_tree,along = 3)
  
  # # compute the original oob errors for each sample and tree
  # preds_y <- predict(object = rf,
  #                    data = x,
  #                    predict.all = T,
  #                    num.threads=num.threads)$predictions
  # oob_error <- (preds_y-y)^2
  # 
  # node_error_oob <- lapply(1:length(y),FUN=function(i){
  #   node_pred_i <- node_pred
  #   node_pred_i[,!oob[i,]] <- 0
  #   return((node_pred_i-y[i])^2)
  # })
  # node_error_oob <- Reduce("+",node_error_oob)/(mean(oob)*length(y))
  # 
  # #for each tree and sample find the oob error and make it in a matrix with 
  # #the original dimensions
  # oob_real <- matrix(node_error_oob[cbind(as.vector(preds_real),
  #                                         rep(1:rf$num.trees,each=length(y)))],
  #                    nrow=length(y))
  
  # for each marker compute which node a sample would end up in for every 
  # allele and aggregate in a weighted manner
  
  imps <- mclapply(1:ncol(x),mc.cores = num.threads,FUN = function(m){
    #m <<- m
    use_trees <- m_use[m,]
    if(sum(use_trees)==0){
      return(rep(0,nrow(x)))
    }else{
      rf_sub <- suppressWarnings(deforest(rf,!use_trees))
      afreq_m <- af_by_tree[m,,]
      x_m <- x
      
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
        
        y_a <- matrix(node_pred[cbind(as.vector(nodes_a),
                                      rep(which(use_trees),each=length(y)))],
                      nrow=length(y))
        
        error_a <- (y_a-y)^2
        oob_m <- oob[,use_trees,drop=F]
        empty_in_node <- matrix(F,nrow = nrow(x),ncol = max_node)
        for(nt in 1:rf_sub$num.trees){
          
          hit_mat <- cbind(1:nrow(x),nodes_a[,nt])
          hit_mat <- hit_mat[oob_m[,nt],]
          
          error_vec <- error_a[,nt]
          node_error_means <- empty_out
          for(j in 1:nrow(hit_mat)){
            node_error_means[hit_mat[j,2]] <- node_error_means[hit_mat[j,2]]+error_vec[hit_mat[j,1]]
          }
          node_counts <- tabulate(hit_mat[,2],max_node)
          node_error_means <- node_error_means/node_counts
          
          error_a[,nt] <- error_a[,nt]-node_error_means[nodes_a[,nt]]
        }
        
        return(error_a)
      })
      
      a_preds <- abind(a_preds,along = 3)
      
      real_error <- matrix(NA,nrow = nrow(x),ncol = rf_sub$num.trees)
      blended_error <- matrix(0,nrow = nrow(x),ncol = rf_sub$num.trees)
      for(i in 1:length(global_alleles)){
        use_samples <- x[,m]==global_alleles[i]
        real_error[use_samples,] <- a_preds[use_samples,,i]
        
        blended_error <- blended_error+t(t(a_preds[,,i])*afreq_m[,use_trees,drop=F][i,])
      }
      diff_error <- blended_error-real_error
      diff_error[!oob[,use_trees,drop=F]] <- 0
      
      return(rowMeans(diff_error)*mean(m_use_orig[m,]))
    }
  })
  imps <- do.call("rbind",imps)
  
  # return the result matrices
  return(imps)
}

