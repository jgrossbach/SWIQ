# swiq with the per-marker computation in C++ (swiq_importances.cpp).
# Compile the C++ part once per R session before calling swiq():
#   Rcpp::sourceCpp("swiq_importances.cpp")

swiq <- function(y,
                 x,
                 num.trees=10000,
                 num.threads=1,
                 min.bucket=15,
                 imp.blocks=1,
                 max.trees=Inf,
                 min.trees=0,
                 pop.str=NULL,
                 pick.threshold=1,
                 mtry=floor(sqrt(ncol(x)))){
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
                    mtry=mtry,
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
  p <- pbinom(n_pos-1,size = n,prob = 0.5,lower.tail = F)
  return(p)
}

get_importances <- function(x,
                            y,
                            num.trees.block,
                            num.threads,
                            min.bucket,
                            mtry,
                            max.trees,
                            min.trees,
                            pop.str,
                            pick.threshold){
  #impute NAs, including pop structure predictors
  nx <- ncol(x)
  for(i in 1:nx){
    if(any(is.na(x[,i]))){
      x[is.na(x[,i]),i] <- sample(x[!is.na(x[,i]),i],
                                  size=sum(is.na(x[,i])),
                                  replace=T)
    }
  }
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
               mtry=mtry,
               keep.inbag = T,
               always.split.variables = pop.str.names,
               respect.unordered.factors = "partition")
  
  #compute in which trees a sample is oob
  # this has to be computed before deforestation
  oob <- do.call("cbind",rf$inbag.counts)==0
  
  # which markers were used where?
  # same as (1:ncol(x))%in%(treeInfo(rf,nt)$splitvarID+1) for every tree
  m_use_orig <- swiq_split_var_use(rf$forest$child.nodeIDs,
                                   rf$forest$split.varIDs,
                                   ncol(x))
  # Do not consider population structure predictors
  # for which trees are interesting
  m_use_orig[pop.str.ind,] <- F
  
  n_picked <- floor(mean(colSums(m_use_orig[1:nx,])))
  p_pick <- n_picked/ncol(x)
  
  n_picked <- rowSums(m_use_orig[1:nx,])
  med_picked <- median(n_picked)
  
  p_picked <- pbinom(q = n_picked-1,size = rf$num.trees,prob = med_picked/rf$num.trees,lower.tail = F)
  picked <- p_picked<=pick.threshold
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
  storage.mode(preds_real) <- "integer"
  
  # data and alleles in the numeric coding ranger uses internally
  X <- encode_like_ranger(x, rf$forest)
  alleles <- lapply(1:ncol(x), function(m){
    if(is.factor(x[,m])){
      global_alleles <- factor(levels(x[,m]), levels = levels(x[,m]))
    }else{
      global_alleles <- unique(x[,m])
    }
    encode_column(global_alleles, rf$forest$covariate.levels[[m]])
  })
  scale_m <- vapply(1:ncol(x), function(m) mean(m_use_orig[m,]), numeric(1))
  
  # per-marker importances in C++
  imps <- swiq_marker_importances(child_nodeIDs = rf$forest$child.nodeIDs,
                                  split_varIDs = rf$forest$split.varIDs,
                                  split_values = rf$forest$split.values,
                                  is_ordered = rf$forest$is.ordered,
                                  X = X,
                                  y = y,
                                  preds_real = preds_real,
                                  oob = oob,
                                  m_use = m_use,
                                  alleles = alleles,
                                  scale = scale_m,
                                  num_threads = num.threads)
  imps <- imps[1:nx,]
  # return the result matrices
  return(imps)
}

# numeric coding of one column, as predict.ranger() does it
# (characters become factors, factors are recoded to the training levels,
# then data.matrix() turns factors into their integer codes)
encode_column <- function(xx, levels_train){
  if(is.character(xx)){
    xx <- factor(xx)
  }
  if(!is.null(levels_train)){
    new.levels <- setdiff(levels(xx), levels_train)
    xx <- factor(xx, levels = c(levels_train, new.levels), exclude = NULL)
  }
  if(is.factor(xx)){
    return(as.numeric(as.integer(xx)))
  }
  as.numeric(xx)
}

encode_like_ranger <- function(x, forest){
  if(length(colnames(x)) != length(forest$independent.variable.names) ||
     any(colnames(x) != forest$independent.variable.names)){
    x <- x[, forest$independent.variable.names, drop = FALSE]
  }
  X <- vapply(seq_len(ncol(x)), function(j){
    encode_column(x[[j]], forest$covariate.levels[[j]])
  }, numeric(nrow(x)))
  X <- matrix(X, nrow = nrow(x))
  if(any(is.na(X))){
    stop("Missing values in predictors after encoding")
  }
  X
}
