setwd("/cellfile/cellnet/RF_QTL_package/")
source("scripts/swiq.R")
library(qtl2)
iron <- read_cross2( system.file("extdata", "iron.zip", package="qtl2") )
map <- insert_pseudomarkers(iron$gmap, step=1)
pr <- calc_genoprob(iron, map, error_prob=0.002,cores=1)
Xcovar <- get_x_covar(iron)
out <- scan1(pr, iron$pheno, Xcovar=Xcovar)
out_nocov <- scan1(pr, iron$pheno, Xcovar=Xcovar)
g <- do.call("cbind",iron$geno)
#g[g==0] <- NA
par(mar=c(5.1, 4.1, 1.1, 1.1))
ymx <- maxlod(out) # overall maximum LOD score
plot(out, map, lodcolumn=1, col="slateblue", ylim=c(0, ymx*1.02))
plot(out, map, lodcolumn=2, col="violetred", add=TRUE)
legend("topleft", lwd=2, col=c("slateblue", "violetred"), colnames(out), bg="gray90")

operm <- scan1perm(pr, iron$pheno, Xcovar=Xcovar, n_perm=5000,cores=10)

summary(operm, alpha=c(0.2, 0.05,0.01))
which(out[,1]>=3.33)

p <- swiq(y = iron$pheno[,1],
          x=g,
          num.trees = 20000,
          num.threads = 5,
          imp.blocks = 4,
          pop.str = NULL,min.bucket = 15)
#p2 <- swiq(y = iron$pheno[,2],x=g,num.trees = 100000,num.threads = 5)

# g_fac <- as.data.frame(g)
# for(i in 1:ncol(g_fac)){
#   g_fac[,i] <- as.factor(g_fac[,i])
# }
# p3 <- swiq(y = iron$pheno[,1],x=g_fac,num.trees = 10000,num.threads = 10,imp.blocks =10,pop.str = Xcovar)

hits <- colnames(g)[which(p<0.05)]
plot(y=out[,1],x=1:nrow(out),type = "h",col=rgb(rownames(out)%in%hits,0,rownames(out)%in%setdiff(colnames(g),hits),1))
abline(h=3.52)

res_mat <- cbind(out[colnames(g),1],p)

which(res_mat[,1]>=3.32&p>0.05)

p_aov <- apply(g,2,FUN=function(x){
  aov_res <- summary(aov(iron$pheno[,1]~as.factor(x)))
  aov_res <- aov_res[[1]][1,5]
  return(aov_res)
})

y <- iron$pheno[,1]
x <- g
boxplot(y~g[,6])
