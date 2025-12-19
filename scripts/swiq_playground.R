set.seed(15)
y <- rnorm(1000)
afreq <- rep(seq(100,500,100),each=40)
x <- matrix(0,ncol=length(afreq),nrow=1000)
for(i in 1:ncol(x)){
  x[sample(1:nrow(x),size = afreq[i]),i] <- 1
}
# y[x[,20]==0] <- y[x[,20]==0]+1
# y[x[,41]==0] <- y[x[,41]==0]-1
# y[x[,81]==0] <- y[x[,81]==0]+1
# y[x[,199]==0] <- y[x[,199]==0]-1
x <- as.data.frame(x)
rf <- ranger(y = y,
             x = x,
             importance="permutation",
             local.importance = T,
             num.trees = 10000,
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
n_pos <- rowSums(tst>0,na.rm=T)
p3 <- pbinom(n_pos-1,size = n,prob = 0.5,lower.tail = F)

p4 <- swiq(y = y,x = x,num.trees = 20000,num.threads = 0,min.bucket = 15,correct_alleles = F,rm.overtraining = T,NA.correction = F,min_use = 2000)
p5 <- swiq(y = y,x = x,num.trees = 20000,num.threads = 0,min.bucket = 15,correct_alleles = F,rm.overtraining = T,NA.correction = F)
p6 <- swiq(y = y,x = x,num.trees = 20000,num.threads = 0,min.bucket = 15,correct_alleles = F,rm.overtraining = F)

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
