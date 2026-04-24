setwd("/cellfile/cellnet/RF_QTL_package/")
library(RandomForestExtended)
library(RFQTL)
load("data/eQtlMappingData160728.RData")
load("data/eQTL_results160831.RData")
g0 <- genotype
g0[is.na(genotype)] <- 0
#i <- 2820
pmat <- matrix(NA,nrow = nrow(phenotype),ncol = 3593)


startPoint <- which(apply(is.na(pmat),1,any))[1]
for(i in startPoint:nrow(phenotype)){
  y <- phenotype[i,]
  use_idx <- which(!is.na(y))
  y <- y[use_idx]
  x <- g0[use_idx,]
  #y[g0[,1140]==1] <- y[g0[,1140]==1]+1
  rf <- RandomForestExtended::randomForest(y=y,x=x,importance=TRUE,keep.inbag = T,ntree=30000,localImp = T,do.trace = 1000,nthreads=10)
  # pred_res_all <- abs(rf$y - rf$predicted)
  # pred_res_m <- matrix(NA,nrow=3593,ncol=length(y))
  # for(m in 1:nrow(pred_res_m)){
  #   x_mut <- x
  #   x_mut[,m] <- 1-x_mut[,m]
  #   pred_ind <- predict(object = rf,newdata = x_mut,predict.all = T)$individual
  #   pred_ind[rf$inbag>0] <- NA
  #   pred_agg <- rowMeans(pred_ind,na.rm = T)
  #   pred_res <- abs(y - pred_agg)
  #   pred_res_m[m,] <- pred_res
  # }
  # p_vec <- apply(pred_res_m,1,FUN=function(x){
  #   x <- x-pred_res_all
  #   x[x==0] <- NA
  #   return(wilcox.test(x=x,mu=0,alternative = "greater")$p.value)
  # })
  
  p_vec_binom <- apply(rf$localImportance,1,FUN=function(x){
    n <- length(x)
    n_pos <- sum(x>0)
    n_neg <- n-n_pos
    p <- binom.test(x = c(n_pos,n_neg),n = n,p = 0.5,alternative = "greater")$p.value
    return(p)
  })
  pmat[i,] <- p_vec_binom[1:3593]
  print(i)
}

qmat <- pmat[1:(i-1),]
qmat[] <- p.adjust(qmat,method="fdr")


perm_mat <- matrix(NA,nrow = 10,ncol = 3593)
for(i in 1:nrow(perm_mat)){
  y <- phenotype[i,]
  set.seed(i)
  #y <- rnorm(n = nrow(g0))
  use_idx <- which(!is.na(y))
  y <- y[use_idx]
  x <- g0[use_idx,]
  perm_vec <- sample(x = 1:length(y))
  y <- y[perm_vec]
  x[,3594:3600] <- x[perm_vec,3594:3600]
  #y[g0[,1140]==1] <- y[g0[,1140]==1]+1
  rf_collect <- list()
  marker_used <- matrix(0,nrow = 3593,ncol = length(y))
  marker_used[] <- 0
  for(j in 1:100){
    rf <- RandomForestExtended::randomForest(y=y,x=x,importance=TRUE,keep.inbag = T,ntree=100,localImp = T,do.trace = 1000,nthreads=1)
    used_here <- rf$localImportance!=0
    marker_used[used_here] <- marker_used[used_here]+1
    rf_collect[[j]] <- rf
  }
  big_rf <- do.call("combine",rf_collect)
  var_mat <- big_rf$localImportance
  var_mat[marker_used<(0.4*length(rf_collect))] <- 0
  # pred_res_all <- abs(rf$y - rf$predicted)
  # pred_res_m <- matrix(NA,nrow=3593,ncol=length(y))
  # for(m in 1:nrow(pred_res_m)){
  #   x_mut <- x
  #   x_mut[,m] <- 1-x_mut[,m]
  #   pred_ind <- predict(object = rf,newdata = x_mut,predict.all = T)$individual
  #   pred_ind[rf$inbag>0] <- NA
  #   pred_agg <- rowMeans(pred_ind,na.rm = T)
  #   pred_res <- abs(y - pred_agg)
  #   pred_res_m[m,] <- pred_res
  # }
  # p_vec <- apply(pred_res_m,1,FUN=function(x){
  #   x <- x-pred_res_all
  #   x[x==0] <- NA
  #   return(wilcox.test(x=x,mu=0,alternative = "greater")$p.value)
  # })
  
  p_vec_binom <- apply(var_mat,1,FUN=function(x){
    n <- length(x)
    n_pos <- sum(x>0)
    n_neg <- n-n_pos
    p <- binom.test(x = c(n_pos,n_neg),n = n,p = 0.5,alternative = "greater")$p.value
    return(p)
  })
  perm_mat[i,] <- p_vec_binom[1:3593]
  print(i)
}

perm_q <- perm_mat
perm_q[] <- p.adjust(perm_mat)


barplot(-log10(pmat[1,]))

rf2 <- randomForest::randomForest(y=y,x=x,importance=TRUE,keep.inbag = T,ntree=30000,localImp = T,do.trace = 1000)
rf3 <- randomForest::randomForest(y=y,x=x[,1:3593],importance=TRUE,keep.inbag = T,ntree=30000,localImp = T,do.trace = 1000)


x <- matrix(NA,ncol=10,nrow=112)
x2 <- matrix(NA,ncol=10,nrow=112)
x[,1] <- sample(x = c(0,1),size = nrow(x),replace = T)
for(i in 2:ncol(x)){
  x[,i] <- x[,i-1]
  x[sample(1:nrow(x),size = 1),i] <- sample(c(0,1),size = 1,replace = T)
}
#y[x[,3]==1] <- y[x[,3]==1]+2
x2[,1] <- sample(x = c(0,1),size = nrow(x),replace = T)
for(i in 2:ncol(x2)){
  x2[,i] <- x2[,i-1]
  x2[sample(1:nrow(x),size = 1),i] <- sample(c(0,1),size = 1,replace = T)
}
x <- cbind(x,x2)
perm_mat2 <- matrix(NA,nrow = 2,ncol = 3593)
for(i in 1:2){
  set.seed(i)
  y <- rnorm(n = nrow(x))
  rfx <- RandomForestExtended::randomForest(y=y,x=x,importance=TRUE,keep.inbag = T,ntree=10000,localImp = T,do.trace = 1000)
  p_vec_binom <- apply(rfx$localImportance,1,FUN=function(x){
    n <- length(x)
    n_pos <- sum(x>0)
    n_neg <- n-n_pos
    p <- binom.test(x = c(n_pos,n_neg),n = n,p = 0.5,alternative = "greater")$p.value
    return(p)
  })
  perm_mat2[i,] <- p_vec_binom[1:3593]
  print(i)
}

qm <- perm_mat2
qm[] <- p.adjust(perm_mat2,method="fdr")

local_imp <- rfx$localImportance
local_imp[] <- 0
n <- ncol(local_imp)
#r2_all <- (rfx$predicted-rfx$y)^2
pred_all <- predict(object = rfx,newdata = x,predict.all = T)$individual
pred_all[rfx$inbag!=0] <- NA
r2_all <- (pred_all-rfx$y)^2
for(m in 1:nrow(local_imp)){
  for(i in 1:10){
    x_m <- x
    x_m[,m] <- sample(x_m[,m],size = n)
    preds <- predict(object = rfx,newdata = x_m,predict.all = T)$individual
    preds[rfx$inbag!=0] <- NA
    r2 <- (preds-y)^2
    r2_delta <- r2-r2_all
    local_imp[m,] <- rowMeans(r2_delta,na.rm=T)
  }
}

local_imp <- rfx$localImportance
local_imp[] <- 0
n <- ncol(local_imp)
pred_all <- predict(object = rfx,newdata = x,predict.all = T)$individual
pred_all[rfx$inbag!=0] <- NA
pred_all <- rowMeans(pred_all,na.rm=T)
r2_all <- (rfx$y-pred_all)^2
for(m in 1:nrow(local_imp)){
  for(i in 1:10){
    x_m <- x
    x_m[,m] <- sample(x_m[,m],size = n)
    preds <- predict(object = rfx,newdata = x_m,predict.all = T)$individual
    preds[rfx$inbag!=0] <- NA
    preds <- rowMeans(preds,na.rm = T)
    r2 <- (preds-rfx$y)^2
    local_imp[m,] <- local_imp[m,]+(r2-r2_all)
  }
}


sv <- perm_mat2
avg_pop <- sapply(1:3593,FUN=function(i){
  mean(rfx$forest$nodepops[rfx$forest$bestvar==i])
})

oob_and_marker <- sapply(1:112,FUN=function(i){
  oob <- rfx$inbag[i,]==0
  marker <- colSums(rfx$forest$bestvar==3363)>0
  return(sum(oob&marker))
})

marker_used <- rfx$localImportance
marker_used[] <- 0
l_collect <- list()
for(i in 1:100){
  rfx <- RandomForestExtended::randomForest(y=y,x=x,importance=TRUE,keep.inbag = T,ntree=10000,localImp = T,do.trace = 1000,nthreads=10)
  used_here <- rfx$localImportance!=0
  marker_used[used_here] <- marker_used[used_here]+1
  l_collect[[i]] <- rfx
  print(i)
}
big_rf <- do.call("combine",l_collect)

p_1 <- apply(big_rf$localImportance,1,FUN=function(x){
  n <- length(x)
  n_pos <- sum(x>0)
  n_neg <- n-n_pos
  p <- binom.test(x = c(n_pos,n_neg),n = n,p = 0.5,alternative = "greater")$p.value
  return(p)
})

p_2 <- apply(rfx$localImportance,1,FUN=function(x){
  n <- length(x)
  n_pos <- sum(x>0)
  n_neg <- n-n_pos
  p <- binom.test(x = c(n_pos,n_neg),n = n,p = 0.5,alternative = "greater")$p.value
  return(p)
})

p_3 <- apply(rfx$localImportance,1,FUN=function(x){
  n <- length(x)
  n_pos <- sum(x>0)
  n_neg <- n-n_pos
  p <- binom.test(x = c(n_pos,n_neg),n = n,p = 0.5,alternative = "greater")$p.value
  return(p)
})

var_mat <- big_rf$localImportance
var_mat[marker_used<20] <- 0
p_4 <- apply(var_mat,1,FUN=function(x){
  n <- length(x)
  n_pos <- sum(x>0)
  n_neg <- n-n_pos
  p <- binom.test(x = c(n_pos,n_neg),n = n,p = 0.5,alternative = "greater")$p.value
  return(p)
})

qtl_rf <- function(y,x,ntree,inclusion,blocksize=100){
  use_idx <- which(!is.na(y))
  y <- y[use_idx]
  x <- x[use_idx,]
  n_blocks <- ceiling(ntree/blocksize)
  rf_list <- list()
  marker_used <- matrix(0,nrow = ncol(x),ncol = length(y))
  marker_used[] <- 0
  for(i in 1:n_blocks){
    x_i <- apply(x,2,FUN=function(x){
      na_bol <- is.na(x)
      if(any(na_bol)){
        x[na_bol] <- sample(x = x[!na_bol],size = sum(na_bol),replace = T)
      }
      return(x)
    })
    rf <- RandomForestExtended::randomForest(y=y,x=x_i,importance=TRUE,keep.inbag = T,ntree=blocksize,localImp = T,nthreads=1,nPerm = 10)
    marker_used[rf$localImportance!=0] <- marker_used[rf$localImportance!=0]+1
    rf_list[[i]] <- rf
  }
  
  big_rf <- do.call("combine",rf_list)
  var_mat <- big_rf$localImportance
  var_mat[marker_used<(inclusion*length(rf_list))] <- 0
  p_vec <- apply(var_mat,1,FUN=function(x){
    n <- length(x)
    n_pos <- sum(x>0)
    n_neg <- n-n_pos
    p <- binom.test(x = c(n_pos,n_neg),n = n,p = 0.5,alternative = "greater")$p.value
    return(p)
  })
  return(p_vec)
}

y_perm <- phenotype[100,]
y_perm <- y_perm[sample(1:length(y_perm))]
p1 <- qtl_rf(y = phenotype[100,],x = genotype[,1:3593],ntree = 20000,inclusion = 0.3)
p2 <- qtl_rf(y = phenotype[100,],x = genotype[,1:3593],ntree = 20000,inclusion = 0)
p3 <- qtl_rf(y = phenotype[100,],x = genotype[,1:3593],ntree = 50000,inclusion = 0.3)
p4 <- qtl_rf(y = phenotype[100,],x = genotype[,1:3593],ntree = 50000,inclusion = 0)
p5 <- qtl_rf(y = phenotype[100,],x = genotype[,1:3593],ntree = 100000,inclusion = 0)
p6 <- qtl_rf(y = y_perm,x = genotype[,1:3593],ntree = 50000,inclusion = 0)
p7 <- qtl_rf2(y = y_perm,x = genotype[,1:3593],ntree = 50000)
p8 <- qtl_rf2(y = phenotype[100,],x = genotype[,1:3593],ntree = 10000)
p9 <- qtl_rf2(y = phenotype[100,],x = genotype[,1:3593],ntree = 1000,n_geno_perm = 1)
p10 <- qtl_rf2(y = y_perm,x = genotype[,1:3593],ntree = 10000)
p11 <- qtl_rf2(y = y_perm,x = genotype[,1:3593],ntree = 10000)

qtl_rf2 <- function(y,x,ntree,blocksize=100,n_geno_perm=3){
  use_idx <- which(!is.na(y))
  y <- y[use_idx]
  x <- x[use_idx,]
  n_blocks <- ceiling(ntree/blocksize)
  rf_list <- list()
  for(i in 1:n_blocks){
    x_i <- apply(x,2,FUN=function(x){
      na_bol <- is.na(x)
      if(any(na_bol)){
        x[na_bol] <- sample(x = x[!na_bol],size = sum(na_bol),replace = T)
      }
      return(x)
    })
    rf <- RandomForestExtended::randomForest(y=y,x=x_i,importance=F,keep.inbag = T,ntree=blocksize,localImp = T,nPerm = 1)
    rf_list[[i]] <- rf
  }
  
  big_rf <- do.call("combine",rf_list)
  
  # local_imp <- matrix(0,nrow = ncol(x),ncol = length(y))
  # for(i in 1:n_geno_perm){
  #   x_i <- apply(x,2,FUN=function(x){
  #     na_bol <- is.na(x)
  #     if(any(na_bol)){
  #       x[na_bol] <- sample(x = x[!na_bol],size = sum(na_bol),replace = T)
  #     }
  #     return(x)
  #   })
  #   for(m in 1:ncol(x)){
  #     x_m <- x_i
  #     x_m[,m] <- 1-x_m[,m]
  #     preds <- predict(object = big_rf,newdata = x_m,predict.all = T)$individual
  #     preds[big_rf$inbag!=0] <- NA
  #     ind_preds <- rowMeans(preds,na.rm = T)
  #     perf_diff <- abs(ind_preds-big_rf$y)-abs(big_rf$predicted-big_rf$y)
  #     local_imp[m,] <- perf_diff+local_imp[m,]
  #   }
  # }
  # local_imp <- local_imp/n_geno_perm
  local_imp <- scale(big_rf$localImportance,center = T,scale=F)
  p_vec <- apply(local_imp,1,FUN=function(x){
    n <- length(x)
    n_pos <- sum(x>0)
    n_neg <- n-n_pos
    p <- binom.test(x = c(n_pos,n_neg),n = n,p = 0.5,alternative = "greater")$p.value
    return(p)
  })
  return(p_vec)
}

local_pred <- matrix(NA,nrow=100,ncol=length(y))
for(i in 1:nrow(local_pred)){
  x_i <- apply(x,2,FUN=function(x){
    na_bol <- is.na(x)
    if(any(na_bol)){
      x[na_bol] <- sample(x = x[!na_bol],size = sum(na_bol),replace = T)
    }
    return(x)
  })
  x_i[,470] <- x_i[sample(1:nrow(x),size=nrow(x)),470]
  preds <- predict(object = big_rf,newdata = x_i,predict.all = T)$individual
  preds[big_rf$inbag!=0] <- NA
  ind_preds <- rowMeans(preds,na.rm = T)
  local_pred[i,] <- ind_preds
}
combined_pred <- colMeans(local_pred)
perf_diff <- abs(combined_pred-big_rf$y)-abs(big_rf$predicted-big_rf$y)
hist(perf_diff)

local_pred2 <- matrix(NA,nrow=100,ncol=length(y))
for(i in 1:nrow(local_pred2)){
  x_i <- apply(x,2,FUN=function(x){
    na_bol <- is.na(x)
    if(any(na_bol)){
      x[na_bol] <- sample(x = x[!na_bol],size = sum(na_bol),replace = T)
    }
    return(x)
  })
  x_i[,1114] <- 1-x_i[,1114]
  preds <- predict(object = big_rf,newdata = x_i,predict.all = T)$individual
  preds[big_rf$inbag!=0] <- NA
  ind_preds <- rowMeans(preds,na.rm = T)
  local_pred2[i,] <- ind_preds
}
combined_pred2 <- colMeans(local_pred2)
perf_diff2 <- abs(combined_pred2-big_rf$y)-abs(big_rf$predicted-big_rf$y)


pmat <- matrix(NA,nrow = nrow(phenotype),ncol = ncol(x))
for(i in 1:nrow(pmat)){
  y <- phenotype[i,]
  #y <- y[sample(1:length(y))]
  rf <- ranger(y = y, x = data.frame(genotype[,1:3593]),importance="permutation",local.importance = T,num.trees = 10000,num.threads = 10,min.bucket = 15)
  p_vec <- apply(t(rf$variable.importance.local),1,FUN=function(x){
    n <- length(x)
    n_pos <- sum(x>0)
    n_neg <- n-n_pos
    p <- binom.test(x = c(n_pos,n_neg),n = n,p = 0.5,alternative = "greater")$p.value
    return(p)
  })
  pmat[i,] <- p_vec
  print(i)
}
qmat <- pmat
qmat[] <- p.adjust(qmat,method="fdr")

perms <- qmat

save(perms,qmat,qv,file = "data/qValues.RData")

setwd("/cellfile/cellnet/phosphoQTL/")
source("lib/pvalues.R")
source("lib/QTL_draft.R")
load("data/plotPack.RData")

rownames(qmat) <- rownames(phenotype)

geno <- read.table("data/genotype_for_mapping.tsv",header=T)
chrVec <- geno[,1]
geno <- geno[,4:ncol(geno)]
geno <- t(geno)
fdrs <- c(1,5,10,15,20,25) 
QtlList <- lapply(fdrs,FUN=function(fdr){
  QTLgrouper(pmat = qmat,sigThreshold = fdr/100,corThreshold = 0.8,distThreshold = 9,genotype = geno,chrVec = chrVec)
})
names(QtlList) <- paste0("FDR",fdrs)

#colors
library(RColorBrewer)
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]
names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")

names(chr_length) <- gsub(pattern = "chr",replacement = "",x = names(chr_length))
png(filename = "graph/eQTL/eQTL_FDR10.png",width = 960,height = 960)
par(cex=1.5,cex.lab=2,mar=c(5,5,4,2)+0.1)
QTLplotter(QTLlist = QtlList$FDR10,targetLocs = geneLocs[rownames(qv)],main="",predictorLocs = markerLocs,chrLen = chr_length,xlab = "eQTL position",ylab = "target-gene position",qtl.lwd=5,labcex = 1,col=adjustcolor(col = colors["eQTL"],alpha.f = 0.8),chrCol = rgb(0,0,0,0.3),chrLty=1)
dev.off()


# get gene prosition (chr coordinates)
gff <- read.table(file = "../phosphoQTL/Saccharomyces_cerevisiae/RM/genes_RM_CDSonly_refined.gtf",stringsAsFactors = F)
CDS_position<-t(sapply(unique(gff[,9]),function(gene){
  sel<-gff[,9]==gene
  pos<-c(gff[sel,4],gff[sel,5])
  pos<-range(pos)
  chr<-gff[sel,1][1]
  return(list(chr,pos[1],pos[2]))
}))
gene<-CDS_position
gene<-data.frame(chr=unlist(gene[,1]),
                 start=unlist(gene[,2]),
                 end=unlist(gene[,3]),
                 stringsAsFactors = F)
gene$med <- rowMeans(cbind(gene$start,gene$end))

# get gene position (relative to genetic markers)
g <- read.table("data/genotype_for_mapping.tsv",header=T,check.names = F)
colnames(g)[1] <- "chr"
geneMarker<-sapply(1:nrow(gene),function(i){
  i <<- i
  before <- max(which(as.character(g$chr)==as.character(gene$chr[i]) &
                        g$end<gene$med[i]))
  if(is.infinite(before)) {before =min(which(as.character(g$chr)==as.character(gene$chr[i])))}
  after <- min(which(as.character(g$chr)==as.character(gene$chr[i]) &
                       genotype$start>gene$med[i]))
  if(is.infinite(after)) {after =max(which(as.character(g$chr)==as.character(gene$chr[i])))}
  return(c(before,after))
})
gene <- cbind(gene,marker_1=geneMarker[1,],marker_2=geneMarker[2,],stringsAsFactors=F)



eQTL <- list(pmat, qmat, QtlList=QtlList)

# identify local QTL
g<-genotype[,4:ncol(genotype)]
g<-t(as.matrix(g))
cor.all <- cor(g,use="pair")

# eQTL
localeQTL<-sapply(eQTL$QtlList$FDR10,function(qtl){
  target <- rownames(phenotype)[qtl$target]
  M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
  if(any(is.infinite(M)))return(FALSE)
  COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
  COR>0.8
})

localeQTLTarget<-unique(sapply(eQTL$QtlList$FDR10[localeQTL],function(x)x$target))
haseQTLTarget <- unique(sapply(eQTL$QtlList$FDR10,function(x)x$target))


NA_p <- t(sapply(1:500,FUN=function(i){
  x <- matrix(NA,ncol=10,nrow=112)
  x[,1] <- sample(x = c(0,1),size = nrow(x),replace = T)
  for(i in 2:ncol(x)){
    x[,i] <- x[,i-1]
    #change all samples for each marker
    x[sample(1:nrow(x),size = 112),i] <- sample(c(0,1),size = 112,replace = T)
  }
  x[1:50,3] <- NA
  y <- rnorm(n = 112)
  p <- swiq(y=y,x=x,min.bucket = 15)
  return(p)
}))

NA_p_corrected <- t(sapply(1:500,FUN=function(i){
  x <- matrix(NA,ncol=10,nrow=112)
  x[,1] <- sample(x = c(0,1),size = nrow(x),replace = T)
  for(i in 2:ncol(x)){
    x[,i] <- x[,i-1]
    #change all samples for each marker
    x[sample(1:nrow(x),size = 112),i] <- sample(c(0,1),size = 112,replace = T)
  }
  x[1:50,3] <- NA
  y <- rnorm(n = 112)
  p <- swiq(y=y,x=x,min.bucket = 15,NA.correction = T)
  return(p)
}))

rf <- ranger(y = phenotype[10,],
             x = data.frame(genotype[,1:3593]),
             importance="permutation",
             local.importance = T,
             num.trees = 10000,
             num.threads = 0,
             min.bucket = 15,
             na.action = "na.learn")


y <- rnorm(n = 112)

NAlist <- extractNAs(genotype[,1:3593])
collect_imps <- matrix(0,nrow = 112,ncol = 3593)
for(i in 1:100){
    #x <- replaceGenoNAs(genotype[,1:3593],NAlist = NAlist)
  x <- genotype[,1:3593]
    rf <- ranger(y = y,
               x = data.frame(x),
               importance="permutation",
               local.importance = T,
               num.trees = 1000,
               num.threads = 0,
               min.bucket = 15,
               na.action = "na.learn")
    collect_imps <- collect_imps+rf$variable.importance.local
}

collect_imps2 <- matrix(0,nrow = 112,ncol = 3593)
for(i in 1:100){
  x <- replaceGenoNAs(genotype = genotype[,1:3593],NAlist = NAlist)
  #x <- genotype[,1:3593]
  rf <- ranger(y = y,
               x = data.frame(x),
               importance="permutation",
               local.importance = T,
               num.trees = 1000,
               num.threads = 0,
               min.bucket = 15,
               na.action = "na.learn")
  collect_imps2 <- collect_imps2+rf$variable.importance.local
}

rf <- ranger(y = y,
             x = data.frame(x),
             importance="permutation",
             local.importance = T,
             num.trees = 10000,
             num.threads = 0,
             min.bucket = 1,
             na.action = "na.learn")

g0 <- genotype
g0[is.na(g0)] <- 0
i <- 0
eQTL_p_perm_imputed <- t(apply(phenotype[1:500,],1,FUN = function(y){
  i <<- i+1
  idx <- sample(1:length(y))
  swiq(y=y[idx],x=g0[,1:3593],min.bucket = 15,num.tree=10000)
}))


i <- 0
eQTL_p_perm_imputed2 <- t(apply(phenotype[1:500,],1,FUN = function(y){
  i <<- i+1
  idx <- sample(1:length(y))
  swiq(y=y[idx],x=genotype[,1:3593],min.bucket = 15,num.tree=10000,multi = T,imp.blocks = 4)
}))

a <- proc.time()
for(i in 1:20000){
  x_i <- apply(x,2,FUN=function(x_m){
    na_bol <- is.na(x_m)
    if(any(na_bol)){
      x_m[na_bol] <- sample(x_m[!na_bol],size = sum(na_bol),replace=T)
    }
    return(x)
  })
  x_i <- data.frame(x_i)
}
b <- proc.time()
b-a

fq <- eQTL_q_noNA
fq[rowSums(qv<0.1)>0,] <- 1

fl <- lapply(fdrs,FUN=function(fdr){
  QTLgrouper(pmat = fq,
             sigThreshold = fdr/100,
             corThreshold = 0.8,
             distThreshold = 9,
             genotype = geno,
             chrVec = chrVec)
})
names(fl) <- paste0("FDR",fdrs)

QTLplotter(QTLlist = fl$FDR10,
           targetLocs = geneLocs[rownames(qv)],
           predictorLocs = markerLocs,
           chrLen = chr_length,
           xlab = "eQTL position",
           ylab = "target-gene position",
           qtl.lwd=5,
           labcex = 1,
           col=adjustcolor(col = colors["eQTL"],alpha.f = 0.8),
           chrCol = rgb(0,0,0,0.3),
           chrLty=1,
           main="swiq results, FDR<10%, only traits with no previous QTL")

gff <- read.table(file = "../phosphoQTL/Saccharomyces_cerevisiae/RM/genes_RM_CDSonly_refined.gtf",stringsAsFactors = F)
CDS_position<-t(sapply(unique(gff[,9]),function(gene){
  sel<-gff[,9]==gene
  pos<-c(gff[sel,4],gff[sel,5])
  pos<-range(pos)
  chr<-gff[sel,1][1]
  return(list(chr,pos[1],pos[2]))
}))
gene<-CDS_position
gene<-data.frame(chr=unlist(gene[,1]),
                 start=unlist(gene[,2]),
                 end=unlist(gene[,3]),
                 stringsAsFactors = F)
gene$med <- rowMeans(cbind(gene$start,gene$end))

g <- read.table("data/genotype_for_mapping.tsv",header=T,check.names = F)
colnames(g)[1] <- "chr"
geneMarker<-sapply(1:nrow(gene),function(i){
  i <<- i
  before <- max(which(as.character(g$chr)==as.character(gene$chr[i]) &
                        g$end<gene$med[i]))
  if(is.infinite(abs(before))) {before =min(which(as.character(g$chr)==as.character(gene$chr[i])))}
  after <- min(which(as.character(g$chr)==as.character(gene$chr[i]) &
                       g$start>gene$med[i]))
  if(is.infinite(after)) {after =max(which(as.character(g$chr)==as.character(gene$chr[i])))}
  return(c(before,after))
})
gene <- cbind(gene,marker_1=geneMarker[1,],marker_2=geneMarker[2,],stringsAsFactors=F)

cor.all <- cor(genotype[,1:3593],use="pair")

# eQTL
localeQTL<-sapply(QtlList$FDR10,function(qtl){
  target <- rownames(phenotype)[qtl$target]
  M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
  if(any(is.infinite(M)))return(FALSE)
  COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
  COR>0.8
})

localeQTLTarget<-unique(sapply(QtlList$FDR10[localeQTL],function(x)x$target))
haseQTLTarget <- unique(sapply(QtlList$FDR10,function(x)x$target))

localeQTL2<-sapply(eQTL_noNA$FDR10,function(qtl){
  target <- rownames(phenotype)[qtl$target]
  M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
  if(any(is.infinite(M)))return(FALSE)
  COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
  COR>0.8
})

localeQTLTarget2 <-unique(unlist(sapply(eQTL_noNA$FDR10[localeQTL2],function(x)x$target)))
haseQTLTarget2 <- unique(sapply(eQTL$FDR10,function(x)x$target))

qtl_p_legacy <- sapply(QtlList$FDR25,FUN=function(qtl){
  pv[qtl$target,qtl$mostSignificantPredictor]
})
top2k_legacy <- QtlList$FDR25[order(qtl_p_legacy)[1:2000]]
localeQTL_legacy<-sapply(top2k_legacy,function(qtl){
  target <- rownames(phenotype)[qtl$target]
  M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
  if(any(is.infinite(M)))return(FALSE)
  COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
  COR>0.8
})
top2k_legacy <- lapply(top2k_legacy,FUN=function(qtl){
  qtl$predictors <- matrix(rep(qtl$mostSignificantPredictor,2),nrow=1,ncol=2)
  colnames(qtl$predictors) <- c("start","end")
  return(qtl)
})

qtl_p_swiq <- sapply(eQTL_noNA$FDR25,FUN=function(qtl){
  eQTL_p_noNA[qtl$target,qtl$mostSignificantPredictor]
})
top2k_swiq <- eQTL_noNA$FDR25[order(qtl_p_swiq)[1:2000]]
localeQTL_swiq<-sapply(top2k_legacy,function(qtl){
  target <- rownames(phenotype)[qtl$target]
  M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
  if(any(is.infinite(M)))return(FALSE)
  COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
  COR>0.8
})
top2k_swiq <- lapply(top2k_swiq,FUN=function(qtl){
  qtl$predictors <- matrix(rep(qtl$mostSignificantPredictor,2),nrow=1,ncol=2)
  colnames(qtl$predictors) <- c("start","end")
  return(qtl)
})

QTLplotter(QTLlist = top2k_legacy,
           targetLocs = geneLocs[rownames(qv)],
           predictorLocs = markerLocs,
           chrLen = chr_length,
           xlab = "eQTL position",
           ylab = "target-gene position",
           qtl.lwd=5,
           labcex = 1,
           col=adjustcolor(col = colors["eQTL"],alpha.f = 0.8),
           chrCol = rgb(0,0,0,0.3),
           chrLty=1,
           main="top 2k qtl, old approach")

QTLplotter(QTLlist = top2k_swiq,
           targetLocs = geneLocs[rownames(qv)],
           predictorLocs = markerLocs,
           chrLen = chr_length,
           xlab = "eQTL position",
           ylab = "target-gene position",
           qtl.lwd=5,
           labcex = 1,
           col=adjustcolor(col = colors["eQTL"],alpha.f = 0.8),
           chrCol = rgb(0,0,0,0.3),
           chrLty=1,
           main="top 2k qtl, swiq")

localeQTL_legacy<-sapply(top2k_legacy,function(qtl){
  target <- rownames(phenotype)[qtl$target]
  M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
  if(any(is.infinite(M)))return(FALSE)
  COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
  COR>0.8
})




localeQTL2<-sapply(eQTL_noNA$FDR10,function(qtl){
  target <- rownames(phenotype)[qtl$target]
  M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
  if(any(is.infinite(M)))return(FALSE)
  COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
  COR>0.8
})

localeQTLTarget2 <-unique(unlist(sapply(eQTL_noNA$FDR10[localeQTL2],function(x)x$target)))
haseQTLTarget2 <- unique(sapply(eQTL$FDR10,function(x)x$target))