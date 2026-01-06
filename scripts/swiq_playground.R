set.seed(15)
y <- rnorm(1000)
afreq <- rep(seq(100,500,100),each=40)
x <- matrix(0,ncol=length(afreq),nrow=1000)
for(i in 1:ncol(x)){
  x[sample(1:nrow(x),size = afreq[i]),i] <- 1
}
y[x[,20]==0] <- y[x[,20]==0]+1
y[x[,41]==0] <- y[x[,41]==0]-1
y[x[,81]==0] <- y[x[,81]==0]+1
y[x[,199]==0] <- y[x[,199]==0]-1
randomEffects <- rnorm(n = 200,sd = 2)
randomEffects[1:200%%10!=0] <- 0
for(m in 1:ncol(x)){
  y[x[,m]==0] <- y[x[,m]==0]+randomEffects[m]
}
x <- as.data.frame(x)
effect_sizes <- sapply(1:ncol(x),FUN=function(m){
  mean(y[x[,m]==0])-mean(y[x[,m]==1])
})

rf <- ranger(y = y,
             x = x,
             importance="permutation",
             local.importance = T,
             num.trees = 5000,
             num.threads = 0,
             min.bucket = 15,keep.inbag = T)
imps <- t(rf$variable.importance.local)
preds_real <- predict(object = rf,data = x,predict.all = T,type="terminalNodes",num.threads=0)$predictions
# m <- 32
# x_m <- x
# x_m[,m] <- 0
# preds_0 <- predict(object = rf,data = x_m,predict.all = T,type="terminalNodes")$predictions
# x_m[,m] <- 1
# preds_1 <- predict(object = rf,data = x_m,predict.all = T,type="terminalNodes")$predictions

# nt <- 5
# tst <- cbind(preds_real[,nt],preds_0[,nt],preds_1[,nt])
# tst[tst[,2]!=tst[,3],]

preds_y <- predict(object = rf,data = x,predict.all = T,num.threads=0)$predictions
oob <- do.call("cbind",rf$inbag.counts)==0
oob_error <- (preds_y-y)^2
oob_error[!oob] <- NA
max_node <- max(preds_real)
node_error_oob <- sapply(1:rf$num.trees,FUN=function(nt){
  nt <<- nt
  out <- rep(NA,max_node)
  preds_real_nt <- preds_real[,nt]
  preds_real_nt[!oob[,nt]] <- NA
  occ_nodes <- unique(sort(preds_real_nt))
  in_node <- outer(preds_real_nt,occ_nodes,"==")
  #in_node[is.na(in_node)] <- F
  errors <- colSums(in_node*oob_error[,nt],na.rm=T)/colSums(in_node,na.rm=T)
  out[occ_nodes] <- errors
  return(out)
})

oob_real <- matrix(oob_error[cbind(as.vector(preds_real),rep(1:rf$num.trees,each=length(y)))],nrow=length(y))
#oob_0 <- matrix(oob_error[cbind(as.vector(preds_0),rep(1:rf$num.trees,each=length(y)))],nrow=length(y))
#oob_1 <- matrix(oob_error[cbind(as.vector(preds_1),rep(1:rf$num.trees,each=length(y)))],nrow=length(y))

oob_real[!oob] <- NA
#oob_0[!oob] <- NA
#oob_1[!oob] <- NA

#exp_err <- (oob_0*sum(x[,m]==0)+oob_1*sum(x[,m]==1))/length(y)-oob_real
#exp_err <- rowMeans(exp_err,na.rm=T)

cor_mat <- matrix(NA,nrow=ncol(x),ncol=nrow(x))
for(m in 1:ncol(x)){
  x_m <- x
  x_m[,m] <- 0
  preds_0 <- predict(object = rf,data = x_m,predict.all = T,type="terminalNodes",num.threads=0)$predictions
  x_m[,m] <- 1
  preds_1 <- predict(object = rf,data = x_m,predict.all = T,type="terminalNodes",num.threads=0)$predictions
  oob_0 <- matrix(oob_error[cbind(as.vector(preds_0),rep(1:rf$num.trees,each=length(y)))],nrow=length(y))
  oob_1 <- matrix(oob_error[cbind(as.vector(preds_1),rep(1:rf$num.trees,each=length(y)))],nrow=length(y))
  exp_err <- (oob_0*sum(x[,m]==0)+oob_1*sum(x[,m]==1))/length(y)-oob_real
  exp_err <- rowMeans(exp_err,na.rm=T)
  cor_mat[m,] <- exp_err
}
sv <- sapply(1:200,FUN=function(i){cor(imps[i,],cor_mat[i,])})
boxplot(sv~afreq)
abline(h=0)

imps_cor <- imps+cor_mat
n <- length(y)
n_pos <- rowSums(imps>0,na.rm=T)
p <- pbinom(n_pos-1,size = n,prob = 0.5,lower.tail = F)
n_pos <- rowSums(imps_cor>0,na.rm=T)
p2 <- pbinom(n_pos-1,size = n,prob = 0.5,lower.tail = F)

plot(-log10(p),-log10(p2),pch=19,col=rgb((1:200)%in%c(20,41,81,199),0,0,0.5))

tst <- remove_overtraining(imps = imps,rf = rf,x = x,y = y)
tst2 <- (tst-imps)*(-1)
n_pos <- rowSums(tst>0,na.rm=T)
p3 <- pbinom(n_pos-1,size = n,prob = 0.5,lower.tail = F)

p4 <- swiq(y = y,x = x,num.trees = 1000,num.threads = 0,min.bucket = 15,correct_alleles = F,rm.overtraining = T,NA.correction = F)
p5 <- swiq(y = y,x = x,num.trees = 10000,num.threads = 0,min.bucket = 15,correct_alleles = F,rm.overtraining = T,NA.correction = F)
p6 <- swiq(y = y,x = x,num.trees = 20000,num.threads = 0,min.bucket = 15,correct_alleles = F,rm.overtraining = T)

marker_used <- sapply(1:rf$num.trees,FUN=function(nt){
  1:200%in%(treeInfo(rf,nt)[,"splitvarID"]+1)
})
ib <- do.call("cbind",rf$inbag.counts)
afreq_bytree <- apply(ib,2,FUN=function(s){
  ind <- rep(1:nrow(x),s)
  colSums((x[ind,])==0)/sum(s)
})
px <- sapply(1:nrow(afreq_bytree),FUN=function(m){
  wilcox.test(afreq_bytree[m,]~marker_used[m,],alternative="greater")$p.value
})
use_mat <- plot(a)

use_freq <- sapply(1:200,FUN=function(m){
  mean(afreq_bytree[m,marker_used[m,]])
})

mpred <- ranger(y=as.factor(marker_used[74,]),x=as.data.frame(t(ib)),num.trees = 1000,num.threads = 0,class.weights = c(0.5,0.5))
PIs <- recompute_imps(rf,x = x,y = y)

tst <- remove_overtraining(imps = imps,rf = rf,x = x,y = y)
tst2 <- (tst-imps)*(-1)
tst2[tst2<0] <- 0
cor_PI <- PIs-tst2


corByTree <- t(sapply(1:rf$num.trees,FUN=function(nt){
  c1 <- cor(y[oob[,nt]],preds_y[oob[,nt],nt])
  c2 <- cor(oob_pred[oob[,nt],nt],y[oob[,nt]])
  return(c(c1,c2))
}))

tI[,1:4] <- tI[,1:4]+1

path_mat <- matrix(NA,nrow=ncol(x),ncol=nrow(tI))
log_mat <- matrix(1,nrow=ncol(x),ncol=nrow(tI))
for(i in 1:nrow(tI)){
  if(!tI$terminal[i]){
    m <- tI$splitvarID[i]
    v <- tI$splitval[i]
    path_mat[,tI[i,"leftChild"]] <- path_mat[,i]
    path_mat[m,tI[i,"leftChild"]] <- v
    path_mat[,tI[i,"rightChild"]] <- path_mat[,i]
    path_mat[m,tI[i,"rightChild"]] <- v
    log_mat[m,tI[i,"rightChild"]] <- -1
  }
}
path_mat_pruned <- path_mat[,tI$terminal]
log_mat_pruned <- log_mat[,tI$terminal]
meta_mat <- cbind(1,tI$nodeID+1,tI$prediction)
meta_mat <- meta_mat[tI$terminal,]

tst <- path_mat_pruned-xi
tst <- tst*log_mat_pruned

colSums(tst>0,na.rm=T)

min_mat <- matrix(-Inf,nrow=ncol(x),ncol=nrow(tI))
max_mat <- matrix(Inf,nrow=ncol(x),ncol=nrow(tI))
for(i in 1:nrow(tI)){
  if(!tI$terminal[i]){
    m <- tI$splitvarID[i]
    v <- tI$splitval[i]
    
    min_mat[,tI[i,"leftChild"]] <- min_mat[,i]
    max_mat[,tI[i,"leftChild"]] <- max_mat[,i]
    max_mat[m,tI[i,"leftChild"]] <- v
    
    min_mat[,tI[i,"rightChild"]] <- min_mat[,i]
    max_mat[,tI[i,"rightChild"]] <- max_mat[,i]
    min_mat[m,tI[i,"rightChild"]] <- v
    
  }
}

min_mat_pruned <- min_mat[,tI$terminal]
max_mat_pruned <- max_mat[,tI$terminal]
meta_mat <- cbind(1,tI$nodeID,tI$prediction)
meta_mat <- meta_mat[tI$terminal,]

rf_mats <- lapply(1:rf$num.trees,FUN=function(nt){
  nt <<- nt
  tI <- treeInfo(object = rf,nt)
  tI[,1:4] <- tI[,1:4]+1
  min_mat <- matrix(-Inf,nrow=ncol(x),ncol=nrow(tI))
  max_mat <- matrix(Inf,nrow=ncol(x),ncol=nrow(tI))
  for(i in 1:nrow(tI)){
    if(!tI$terminal[i]){
      m <- tI$splitvarID[i]
      v <- tI$splitval[i]
      
      min_mat[,tI[i,"leftChild"]] <- min_mat[,i]
      max_mat[,tI[i,"leftChild"]] <- max_mat[,i]
      max_mat[m,tI[i,"leftChild"]] <- v
      
      min_mat[,tI[i,"rightChild"]] <- min_mat[,i]
      max_mat[,tI[i,"rightChild"]] <- max_mat[,i]
      min_mat[m,tI[i,"rightChild"]] <- v
      
    }
  }
  min_mat_pruned <- min_mat[,tI$terminal]
  max_mat_pruned <- max_mat[,tI$terminal]
  meta_mat <- cbind(nt,tI$nodeID,tI$prediction)
  meta_mat <- meta_mat[tI$terminal,]
  return(list(min=min_mat_pruned,
              max=max_mat_pruned,
              meta=meta_mat))
})

rf_min <- lapply(rf_mats,FUN=function(l){
  return(l$min)
})
rf_min <- do.call("cbind",rf_min)

rf_max <- lapply(rf_mats,FUN=function(l){
  return(l$max)
})
rf_max <- do.call("cbind",rf_max)

rf_meta <- lapply(rf_mats,FUN=function(l){
  return(l$meta)
})
rf_meta <- do.call("rbind",rf_meta)

which(colSums(!(rf_min<=xi&rf_max>xi))==0)
preds <- rf_meta[which(colSums(!(rf_min<=xi&rf_max>xi))==0),]

pred_mat <- matrix(NA,ncol=rf$num.trees,nrow=nrow(x))

a <- proc.time()
for(i in 1:nrow(x)){
  xi <- unlist(x[i,])
  preds <- rf_meta[which(colSums(!(rf_min<=xi&rf_max>xi))==0),]
  pred_mat[i,] <- preds[,3]
}
b <- proc.time()
b-a

not_use <- is.infinite(rf_min)&is.infinite(rf_max)
rf_min[not_use,] <- NA
rf_max[not_use,] <- NA

a <- proc.time()
for(i in 1:nrow(x)){
  #xi <- unlist(x[i,])
  preds <- rf_min<=xi&rf_max>xi
  #pred_mat[i,] <- preds[,3]
}
b <- proc.time()
b-a



a <- proc.time()
for(i in 1:ncol(rf_min)){
  n <- which(!is.infinite(rf_min[,i])|!is.infinite(rf_max[,i]))
  min_sub <- rf_min[n,i]
  max_sub <- rf_max[n,i]
  use <- min_sub<=x[n,]&max_sub>x[n,]
}
b <- proc.time()
b-a

# X <- as.matrix(x)  # n × p
# p <- ncol(X)
# L <- ncol(rf_min)
# 
# pred_mat <- matrix(NA_real_, nrow = nrow(X), ncol = rf$num.trees)
# 
# for (j in seq_len(L)) {
#   inside <- rowSums(X >= rf_min[, j] & X < rf_max[, j]) == p
#   pred_mat[inside, rf_meta[j, 1]] <- rf_meta[j, 3]
# }

pred_ranger <- predict(object=rf,predict.all = T,data = x)$predictions


min_arr <- array(rf_min, dim = c(nrow(rf_min), ncol(rf_min), nrow(x)))
max_arr <- array(rf_min, dim = c(nrow(rf_min), ncol(rf_min), 1))
x_arr <- array(x, dim = c(nrow(rf_min), 1, nrow(x)))

min_l <- min_arr<=x_arr

pred_mat <- predict_from_bounds_cpp(
  X        = as.matrix(x),
  rf_min   = rf_min,
  rf_max   = rf_max,
  rf_meta  = rf_meta,
  num_trees = rf$num.trees
)

x_arr <- array(data = 1:2000,dim = c(2,10,100))
y_m <- matrix(data = NA,nrow=2,ncol=100)
x_arr+y_m


x_i <- x
a <- proc.time()
tst <- predict(object = rf,data=x_i,predict.all = T)$predictions
b <- proc.time()
b-a
dim(tst)

x_i <- rbind(x,x)
a <- proc.time()
tst <- predict(object = rf,data=x_i,predict.all = T)$predictions
b <- proc.time()
b-a
dim(tst)

x_i <- rbind(x,x,x)
a <- proc.time()
tst <- predict(object = rf,data=x_i,predict.all = T)$predictions
b <- proc.time()
b-a
dim(tst)

x_i <- rbind(x,x,x,x)
a <- proc.time()
tst <- predict(object = rf,data=x_i,predict.all = T)$predictions
b <- proc.time()
b-a
dim(tst)

x_i <- do.call("rbind",lapply(1:400,FUN=function(i){x}))
a <- proc.time()
tst <- predict(object = rf,data=x_i,predict.all = T,num.threads = 30,save.memory=T)$predictions
b <- proc.time()
b-a
dim(tst)


x_i <- do.call("rbind",lapply(1:ncol(x),FUN=function(m){
  x_0 <- x_1 <- x
  x_0[,m] <- 0
  x_1[,m] <- 1
  return(rbind(x_0,x_1))
}))

rf <- ranger(y = y,
             x = x,
             importance="permutation",
             local.importance = T,
             num.trees = 10000,
             num.threads = 30,
             min.bucket = 15)

a <- proc.time()
tst <- predict(object = rf,data=x_i,predict.all = T,num.threads = 30,save.memory=T,)$predictions
b <- proc.time()
b-a
dim(tst)



a <- proc.time()
tst <- swiq(y = y,x = x,num.trees = 1000,num.threads = 30,min.bucket = 15,rm.overtraining = T)
b <- proc.time()
b-a


a <- proc.time()
preds <- rf_min<=xi&rf_max>xi
b <- proc.time()
b-a

usage <- rowSums(!is.infinite(rf_min)|!is.infinite(rf_max))

xi <- unlist(x[i,])
hit <- 1:ncol(rf_min)
a <- proc.time()
for(m in order(usage,decreasing = T)){
  print(length(hit))
  hit <- hit[rf_min[m,hit]<=xi[m]]
  hit <- hit[rf_max[m,hit]>xi[m]]
}
b <- proc.time()
b-a

a <- proc.time()
pred_by_tree <- lapply(rf_mats[1:200],FUN=function(l){
  usage <- rowSums(!is.infinite(l$min)|!is.infinite(l$max))>0
  x_sub <- t(x[,usage])
  lmin <- l$min[usage,]
  lmax <- l$max[usage,]
  out <- rep(NA,ncol(x_sub))
  for(i in 1:ncol(x_sub)){
    out[i] <- which(colSums(!(lmin<=x_sub[,i]&lmax>x_sub[,i]))==0)
  }
  return(out)
})
b <- proc.time()
b-a

library(snow)

cl <- makeCluster(30, type = "SOCK")
clusterExport(cl = cl,list = c("x","rf"))
invisible(clusterEvalQ(cl, library(ranger)))
a <- proc.time()
out <- parLapply(cl = cl,x = 1:200,fun = function(m){
  x[,m] <- 0
  p0 <- predict(rf,data = x,predict.all = T,type = "terminalNodes")$predictions
  x[,m] <- 1
  p1 <- predict(rf,data = x,predict.all = T,type = "terminalNodes")$predictions
  return(list(p0,p1))
})
b <- proc.time()
b-a
stopCluster(cl)

m_used <- !is.infinite(rf_min)|!is.infinite(rf_max[m,])

m <- 10
keep <- which(m_used[m,])
xi <- x
xi[,m] <- 0
rf_min_s <- rf_min[,keep]
min_arr <- array(data = rf_min_s,dim = c(nrow(rf_min),length(keep),nrow(x)))
min_arr <- aperm(min_arr, perm = c(1, 3, 2))
min_arr <- min_arr-array(data = t(xi),dim = c(ncol(xi),nrow(xi),length(keep)))


norm_pred <- predict(rf,x,predict.all = T,type = "terminalNodes")$predictions
m <- 10


m_used <- !is.infinite(rf_min)|!is.infinite(rf_max[m,])
nodes_0 <- norm_pred
nodes_1 <- norm_pred
keep_m <- which(m_used[m,])
rf_min_m <- rf_min[,keep_m]
out_0 <- matrix(F,nrow = nrow(x),ncol=ncol(m_used))
xm <- x
xm[,m] <- 0
for(i in 1:nrow(x)){
  xi <- unlist(xm[i,])
  keep <- keep_m
  rf_min_s <- rf_min_m
  keep <- keep[colSums(!rf_min_s<=xi)==0]
  rf_max_s <- rf_max[,keep]
  keep <- keep[colSums(!rf_max_s>xi)==0]
  submeta <- rf_meta[keep,]
  nodes_0[i,submeta[,1]] <- submeta[,2]
}
xm[,m] <- 1
for(i in 1:nrow(x)){
  xi <- unlist(xm[i,])
  keep <- keep_m
  rf_min_s <- rf_min_m
  keep <- keep[colSums(!rf_min_s<=xi)==0]
  rf_max_s <- rf_max[,keep]
  keep <- keep[colSums(!rf_max_s>xi)==0]
  submeta <- rf_meta[keep,]
  nodes_1[i,submeta[,1]] <- submeta[,2]
}

m_by_tree <- sapply(1:rf$num.trees,FUN=function(nt){
  (1:ncol(x))%in%(treeInfo(object = rf,tree = nt)$splitvarID+1)
})

use <- m_by_tree[140,]
rf2 <- deforest(object = rf,which.trees = which(!use))
p <- predict(rf,x,predict.all = T)$predictions
p2 <- predict(rf2,x,predict.all = T)$predictions

a <- proc.time()
p <- predict(rf,x,predict.all = T)$predictions
p <- predict(rf,x,predict.all = T)$predictions
b <- proc.time()
b-a

a <- proc.time()
rf2 <- deforest(object = rf,which.trees = which(!use))
p2 <- predict(rf2,x,predict.all = T)$predictions
p2 <- predict(rf2,x,predict.all = T)$predictions
b <- proc.time()
b-a



g <- genotype[,1:3593]
g[is.na(g)] <- 0
g <- as.data.frame(g)

rf <- ranger(y = phenotype[10,],
             x = g,
             importance="permutation",
             local.importance = T,
             num.trees = 500000,
             num.threads = 30,
             min.bucket = 15)
m_by_tree <- sapply(1:rf$num.trees,FUN=function(nt){
  (1:ncol(g))%in%(treeInfo(object = rf,tree = nt)$splitvarID+1)
})

use <- m_by_tree[2566,]

a <- proc.time()
p <- predict(rf,g,predict.all = T)$predictions
p <- predict(rf,g,predict.all = T)$predictions
b <- proc.time()
b-a

a <- proc.time()
rf2 <- deforest(object = rf,which.trees = which(!use))
p2 <- predict(rf2,g,predict.all = T)$predictions
p2 <- predict(rf2,g,predict.all = T)$predictions
b <- proc.time()
b-a
