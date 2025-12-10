#####Swiq#######################################################################
## Description
# swiq: Sample Wise Importances Qtl-mapping
# This function grows a random forest for a trait and then tests each predictor
# for elevated numbers of samples for which the predictor improved the 
# prediction. This is done with a one sided binomial test.
## Parameters:
# y: Dependent variable. If it contains NAs, the samples are removed.
# x: Independent variables. They can contain NAs.
# num.trees: The number of trees that are generated in the RF.
# num.threads: The number of cores that are used to grow the RF. The default 0
#              takes up all available cores
# min.bucket: The number of inbag samples that have to be contained in each 
#             terminal node. When predictors are correlated, small noisy nodes
#             can occur otherwise and lead to false positives.
# min.node.size: How many samples have to be in a node so it can be split?
# NA.correction: Should samples with missing information for a predictor be 
#                excluded for the PI computation for that predictor?
## Output:
# A vector of p-values, one for each predictor.
swiq <- function(y,
                 x,
                 num.trees=10000,
                 num.threads=0,
                 min.bucket=15,
                 min.node.size=NULL,
                 NA.correction=F,
                 multi=F,
                 imp.blocks=4,
                 method="binomial",
                 correction_mode=1,
                 window_size=50,
                 correct_alleles=F,
                 correction_mat=NULL){
  if(!is.data.frame(x)){
    x <- data.frame(x)
  }
  if(!method%in%c("binomial","t-test")){
    stop("method has to be either binomial or t-test!")
  }
  y <- scale(y,center = T,scale = T)[,1]
  if(any(is.na(y))){
    na_bol <- is.na(y)
    y <- y[!na_bol]
    x <- x[!na.bol,]
  }
  imps <- get_imps(x = x,
                   y = y,
                   num.trees = num.trees,
                   num.threads = num.threads,
                   min.bucket = min.bucket,
                   min.node.size = min.node.size,
                   NA.correction = NA.correction,
                   multi = multi,
                   imp.blocks = imp.blocks)
  if(method=="binomial"){
    if(correct_alleles){
      imps[imps==0] <- NA
      imps <- correct_importances(imps = imps,
                                  x = x,
                                  correction_mode=correction_mode,
                                  window_size=window_size,
                                  correction_mat=correction_mat)
    }
    imps <- imps>0
    if(NA.correction){
      imps[is.na(t(x))] <- F
      n <- rowSums(!is.na(t(x)))
    }else{
      n <- ncol(imps)
    }
    
    n_pos <- rowSums(imps,na.rm=T)
    p <- pbinom(n_pos-1,size = n,prob = 0.5,lower.tail = F)
  }
  if(method=="t-test"){
    imps <- t(rf$variable.importance.local)
    if(NA.correction){
      imps[is.na(t(x))] <- NA
    }
    p <- apply(imps,1,FUN=function(x){
      t.test(x,mu=0,alternative="greater")$p.value
    })
    p[is.na(p)] <- 1
  }
  return(p)
}

###QTL_profile##################################################################
#plot the p-values for a single trait
#Arguments:
#           QTLlist: output of QTLgrouper
#           predictorLocs: location of predictors in the genome in absolute 
#                          cummulative positions
#           chrLen: vector containing the length of each chromosome with names
#           xlab: desired label for the x-axis
#           ylab: desired label for the y-axis
#           main: desired main-title
QTL_profile <- function(p_vec,
                        predictorLocs,
                        chrLen,
                        xlab,
                        ylab,
                        main,
                        col="black",
                        labcex=0.6,
                        chrLty=2,
                        chrCol="red"){
  plot.default(y = -log10(p_vec),
               x = predictorLocs,
               main = main,
               xlab = xlab,
               ylab = ylab,
               type = "h",
               col = col,
               xaxt = "n",
               cex.lab=labcex)
  chrLenCum <- sapply(1:length(chrLen),FUN=function(x){
    sum(chrLen[1:x])
  })
  for(i in 1:(length(chrLenCum)-1)){
    abline(v=chrLenCum[i],lty=chrLty,col=chrCol)
  }
  borders <- c(0,chrLenCum)
  for(i in 1:length(chrLen)){
    mtext(names(chrLen[i]),
          side=1,
          at=mean(borders[i:(i+1)]),
          cex=par()$cex*labcex*0.75)
  }
}

#####histogramm with ggplot#####################################################
gg_histogram <- function(x, breaks = 30, color = "steelblue",
                          xlab = "Value", ylab = "Count", title = "Histogram") {
  # Ensure x is numeric
  if (!is.numeric(x)) {
    stop("Input 'x' must be a numeric vector.")
  }
  
  # Create a data frame for ggplot
  df <- data.frame(x = x)
  
  # Build and return the plot
  ggplot(df, aes(x = x)) +
    geom_histogram(bins = breaks, fill = color, color = "black", alpha = 0.8) +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal(base_size = 14)
}

#####barplot with ggplot########################################################
gg_barplot <- function(values, 
                           title = "Barplot", 
                           xlab = "Index", 
                           ylab = "Value", 
                           bar_color = "steelblue") {
  
  # Check if input is numeric
  if(!is.numeric(values)) {
    stop("Input 'values' must be a numeric vector.")
  }
  
  # Create a data frame for ggplot
  df <- data.frame(Index = seq_along(values), Value = values)
  
  # Create barplot using ggplot2
  p <- ggplot(df, aes(x = factor(Index), y = Value)) +
    geom_bar(stat = "identity", fill = bar_color) +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

gg_boxplot <- function(x, y, 
                       title = NULL, 
                       xlab = NULL, 
                       ylab = NULL, 
                       fill = "steelblue",
                       box_labels = NULL,
                       show_outliers = TRUE) {
  # Basic checks
  if (missing(x) | missing(y)) {
    stop("Please provide both x (grouping) and y (numeric vector).")
  }
  if (!is.numeric(y)) {
    stop("y must be numeric.")
  }
  if (!is.factor(x) && !is.character(x)) {
    x <- factor(x)
  }
  
  # Build data frame internally
  df <- data.frame(group = x, value = y)
  
  # Determine y-limits if outliers are hidden
  ylim_vals <- NULL
  if (!show_outliers) {
    # Compute the lower and upper whisker values for each group
    library(dplyr)
    whiskers <- df %>%
      group_by(group) %>%
      summarise(
        lower = quantile(value, 0.25) - 1.5 * IQR(value),
        upper = quantile(value, 0.75) + 1.5 * IQR(value),
        .groups = "drop"
      )
    ylim_vals <- c(min(whiskers$lower), max(whiskers$upper))
  }
  
  # Create plot
  p <- ggplot(df, aes(x = group, y = value)) +
    geom_boxplot(
      fill = fill,
      outlier.shape = if (show_outliers) 19 else NA
    ) +
    theme_minimal()
  
  # Apply adjusted y-limits if needed
  if (!is.null(ylim_vals)) {
    p <- p + coord_cartesian(ylim = ylim_vals)
  }
  
  # Optional labels
  if (!is.null(title)) p <- p + ggtitle(title)
  
  # Use provided box_labels for the x-axis scale breaks
  if (!is.null(box_labels)) {
    if (length(box_labels) != length(levels(x))) {
      stop("The length of 'box_labels' must match the number of groups (levels of x).")
    }
    p <- p + scale_x_discrete(labels = box_labels)
    if (is.null(xlab)) {
      xlab <- "Group"
    }
  }
  
  if (!is.null(xlab))  p <- p + xlab(xlab)
  if (!is.null(ylab))  p <- p + ylab(ylab)
  
  return(p)
}


correct_importances <- function(imps,
                                x,
                                correction_mode,
                                window_size,
                                correction_mat){
  if(correction_mode==1){
    n_0 <- colSums(x==0,na.rm=T)
    n_1 <- colSums(x==1,na.rm=T)
    mean_0 <- rowSums(imps*t(x==0),na.rm=T)/n_0
    mean_1 <- rowSums(imps*t(x==1),na.rm=T)/n_1
    new_mean <- pmin(mean_0,mean_1)
    imps_cor <- (imps-mean_0+new_mean)*t(x==0)+(imps-mean_1+new_mean)*t(x==1)
    return(imps_cor)
  }
  if(correction_mode==2){
    n_0 <- colSums(x==0,na.rm=T)
    n_1 <- colSums(x==1,na.rm=T)
    mean_0 <- rowSums(imps*t(x==0),na.rm=T)/n_0
    mean_1 <- rowSums(imps*t(x==1),na.rm=T)/n_1
    major <- rep(0,ncol(x))
    major[n_1>n_0] <- 1
    major_af <- pmax(n_0,n_1)/(n_0+n_1)
    af_r <- rank(major_af,ties.method = "first")
    af_r <- ncol(x)-af_r+1
    window_edges <- cbind(af_r-ceiling(window_size/2),
                          af_r+ceiling(window_size/2))
    
    too_small <- window_edges[,1]<1
    window_edges[too_small,] <- window_edges[too_small,] + 
      abs(window_edges[too_small,1])+1
    too_big <- window_edges[,2]>ncol(x)
    window_edges[too_big,] <- window_edges[too_big,] - 
      (window_edges[too_big,2]-ncol(x))
    major_means <- cbind(mean_0,mean_1)[cbind(1:ncol(x),major+1)]
    imps_cor <- imps
    for(m in 1:ncol(x)){
      i_m <- imps[m,]
      window_markers <- setdiff(window_edges[m,1]:window_edges[m,2],m)
      correction_factor <- mean(major_means[window_markers],na.rm=T)
      correction_factor <- max(correction_factor,0)
      i_m[which(x[,m]==major[m])] <- i_m[which(x[,m]==major[m])]-correction_factor
      imps_cor[m,] <- i_m
    }
    return(imps_cor)
  }
  if(correction_mode==3){
    imps_cor <- imps
    for(m in 1:ncol(x)){
      a_l <- sort(unique(x[,m]))
      for(a_i in 1:length(a_l)){
        a <- a_l[a_i]
        imps_cor[m,which(x[,m]==a)] <- imps_cor[m,which(x[,m]==a)]-
          correction_mat[m,a_i]
      }
    }
    return(imps_cor)
  }

}



correct_importances_generalized <- function(imps,x,exclude,freq_min){
  x_out <- x
  if(!is.null(exclude)){
    x <- x[,-exclude]
  }
  for(m in 1:ncol(x)){
    x_levels <- sort(unique(x[,m]))
    na_pos <- is.na(x[,m])
    matches <- which(outer(X = x[,m],Y = x_levels,FUN = "=="),arr.ind=T)
    matches <- matches[order(matches[,1]),]
    alleles <- matches[,2]
    x[!na_pos,m] <- alleles
  }
  max_allele <- max(x,na.rm=T)
  n_a <- matrix(NA,nrow = max_allele,ncol = ncol(x))
  
  mean_a <- n_a
  for(a in 1:max_allele){
    n_a[a,] <- colSums(x==a)
    mean_a[a,] <- rowMeans(imps*t(x==a))
  }
  rm_alleles <- which(n_a<freq_min,arr.ind=T)
  for(i in 1:nrow(rm_alleles)){
    m <- rm_alleles[i,2]
    a <- rm_alleles[i,1]
    imps[m,x[,m]==a] <- NA
  }
  mean_a[rm_alleles] <- NA
  new_mean <- mean_a*n_a
  #new_mean <- (mean_0*n_0+mean_1*n_1)/(n_0+n_1)
  new_mean <- apply(mean_a,2,min,na.rm=T)
  clean_imps <- lapply(1:max_allele,FUN=function(a){
    (imps-mean_a[a,])*t(x==a)
  })
  clean_imps <- Reduce(f = "+",clean_imps)
  imps_cor <- clean_imps+new_mean
  imps_cor[is.na(imps)] <- NA
  return(imps_cor)
}

get_imps <- function(x,
                     y,
                     num.trees=10000,
                     num.threads=0,
                     min.bucket=15,
                     min.node.size=NULL,
                     NA.correction=F,
                     multi=F,
                     imp.blocks=4){
  if(multi){
    imps <- matrix(0,nrow = nrow(x),ncol = ncol(x))
    for(i in 1:imp.blocks){
      x_i <- apply(x,2,FUN=function(x_m){
        na_bol <- is.na(x_m)
        if(any(na_bol)){
          x_m[na_bol] <- sample(x_m[!na_bol],size = sum(na_bol),replace=T)
        }
        return(x_m)
      })
      x_i <- data.frame(x_i)
      rf <- ranger(y = y,
                   x = x_i,
                   importance="permutation",
                   local.importance = T,
                   num.trees = ceiling(num.trees/imp.blocks),
                   num.threads = num.threads,
                   min.bucket = min.bucket,
                   min.node.size = min.node.size)
      imps <- imps+rf$variable.importance.local
    }
    imps <- t(imps)
  }else{
    rf <- ranger(y = y,
                 x = x,
                 importance="permutation",
                 local.importance = T,
                 num.trees = num.trees,
                 num.threads = num.threads,
                 min.bucket = min.bucket,
                 min.node.size = min.node.size)
    imps <- t(rf$variable.importance.local)
    
  }
  return(imps)
}

estimate_correction <- function(x,
                                num.trees=10000,
                                num.threads=0,
                                min.bucket=15,
                                min.node.size=NULL,
                                NA.correction=F,
                                multi=F,
                                imp.blocks=4,
                                n_sim=10){
  if(!is.data.frame(x)){
    x <- data.frame(x)
  }
  max_allele <- max(apply(x,2,FUN=function(m){length(unique(m))}))
  allele_bias <- matrix(0,nrow = ncol(x),ncol = max_allele)
  for(i in 1:n_sim){
    y <- rnorm(nrow(x))
    y <- scale(y,center = T,scale = T)[,1]
    set.seed(as.numeric(Sys.time()))
    imps <- get_imps(x = x,
                     y = y,
                     num.trees = num.trees,
                     num.threads = num.threads,
                     min.bucket = min.bucket,
                     min.node.size = min.node.size,
                     NA.correction = NA.correction,
                     multi = multi,
                     imp.blocks = imp.blocks)
    for(a_i in 0:(max_allele-1)){
      allele_bias[,a_i+1] <- allele_bias[,a_i+1]+rowSums((t(x)==a_i)*imps)/rowSums(t(x)==a_i)
    }
    
  }
  allele_bias <- allele_bias/n_sim
  return(allele_bias)
}

named_list <- function(name_vec){
  l <- lapply(name_vec,get)
  names(l) <- name_vec
  return(l)
}

gg_dotplot <- function(x, y,
                      col = "black",
                      main = NULL,
                      pch = 19,
                      xlab = NULL,
                      ylab = NULL,
                      xlim = NULL,
                      ylim = NULL) {
  
  # Load required package
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  
  # Build data frame
  df <- data.frame(x = x, y = y)
  
  # Start plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(color = col, shape = pch)
  
  # Labels
  if (!is.null(main)) p <- p + ggplot2::ggtitle(main)
  if (!is.null(xlab)) p <- p + ggplot2::xlab(xlab)
  if (!is.null(ylab)) p <- p + ggplot2::ylab(ylab)
  
  # Axis limits
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  }
  
  return(p)
}
