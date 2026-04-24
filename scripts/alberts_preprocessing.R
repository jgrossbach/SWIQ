library(ranger)
library(readxl)

traits <- read.table("data/albert2018/SI_Data_01_expressionValues.txt",header = T,quote = "",row.names = 1)
traits <- t(traits)

genotype <- read.table("data/albert2018/SI_Data_03_genotypes.txt",header=T)
genotype[genotype==(-1)] <- 0
gcor <- sapply(1:(ncol(genotype)-1),FUN=function(i){
  cor(genotype[,i],genotype[,i+1])
})

chrvec <- gsub(x = colnames(genotype),pattern = "[.].*",replacement = "")
chrs <- unique(chrvec)

markerGroups <- lapply(chrs,FUN=function(chr){
  chr <<- chr
  chr_snps <- which(chrvec==chr)
  # allCor <- cor(genotype[,chr_snps],use="pair")^2
  # d <- as.dist(1-allCor)
  d <- dist(t(genotype[,chr_snps]),method = "manhattan")
  #d <- as.matrix(d)
  h <- hclust(d)
  genotype2group <- cutree(h,h = 10)
  return(genotype2group)
})

running_n <- 0
for(i in 1:length(markerGroups)){
  markerGroups[[i]] <- markerGroups[[i]]+running_n
  running_n <- max(markerGroups[[i]])
}

markerGroups <- do.call("c",markerGroups)
genotype_red <- sapply(unique(markerGroups),FUN=function(m){
  m <<- m
  apply(genotype[,markerGroups==m,drop=F],1,median)
})
genotype_red[genotype_red==0.5] <- 0

covariates <- read_excel("data/albert2018/SI_Data_02_covariates.xlsx")
qtl_alberts <- read_excel("data/albert2018/SI_Data_04_eQTL.xlsx")
qtl_alberts$pmarker <- gsub(pattern = ":",replacement = ".",x = qtl_alberts$pmarker)
qtl_alberts$pmarker <- gsub(pattern = "/",replacement = ".",x = qtl_alberts$pmarker)

qtl_alberts$CI.l <- gsub(pattern = ":",replacement = ".",x = qtl_alberts$CI.l)
qtl_alberts$CI.l <- gsub(pattern = "/",replacement = ".",x = qtl_alberts$CI.l)

qtl_alberts$CI.r <- gsub(pattern = ":",replacement = ".",x = qtl_alberts$CI.r)
qtl_alberts$CI.r <- gsub(pattern = "/",replacement = ".",x = qtl_alberts$CI.r)


batches <- sort(unique(covariates$batch))
batch_est <- sapply(batches,FUN=function(b){
  rowMeans(traits[,covariates$batch==b,drop=F])
})
colnames(batch_est) <- batches
batch_est <- batch_est-rowMeans(batch_est,na.rm=T)
traits <- traits-batch_est[,covariates$batch]
traits <- traits-rowMeans(traits,na.rm=T)
traits <- traits/apply(traits,1,sd,na.rm=T)

alberts_data <- list(genotype_full=genotype,
                     genotype=genotype_red,
                     marker_groups=markerGroups,
                     phenotype=traits)

saveRDS(object = alberts_data,file = "data/albert_mapping_data.rds")

#test1
g <- "YAL021C"
y <- traits[g,]
x <- genotype_red

fc <- sapply(1:ncol(x),FUN=function(m){
  mean(y[x[,m]==1])-mean(y[x[,m]==0])
})

# y_pred <- ranger(y = y,
#                  x = as.data.frame(covariates[,3]),
#                  importance="none",
#                  local.importance = F,
#                  num.trees = 5000,
#                  num.threads = 15,
#                  min.bucket = 15,
#                  keep.inbag = T)

a <- proc.time()
tst <- mclapply(X = 1:10,mc.cores = 10,FUN=function(i){
  
  p <- swiq(y = y,
            x = x,
            num.trees = 20000,
            num.threads = 1,
            min.bucket = 15,
            imp.blocks = 1,
            max.trees = Inf,
            min.trees = 100,
            pick.threshold = 0.1)
  
})
b <- proc.time()
(b-a)[3]

a <- proc.time()
p <- swiq(y = y,
          x = x,
          num.trees = 20000,
          num.threads = 1,
          min.bucket = 15,
          imp.blocks = 1,
          max.trees = Inf,
          min.trees = 100,
          pick.threshold = 0.1)
b <- proc.time()
(b-a)[3]

q <- p.adjust(p,method = "fdr")
p_mapped <- q[markerGroups]
names(p_mapped) <- colnames(genotype)

ref_hits <- qtl_alberts[qtl_alberts$gene==g,]
ref_hits <- lapply(1:nrow(ref_hits),FUN=function(i){
  m1 <- which(names(p_mapped)==unlist(ref_hits[i,"CI.l"]))
  m2 <- which(names(p_mapped)==unlist(ref_hits[i,"CI.r"]))
  return(m1:m2)
})
#ref_hits <- sort(unlist(ref_hits))
plot(1:length(p_mapped),-log10(p_mapped),type="h")
#points(x = which(names(p_mapped)%in%ref_hits),y=rep(max(-log10(p_mapped)),length(ref_hits)))
lapply(ref_hits,FUN=function(v){
  rect(xleft = v[1],ybottom = 0,xright = tail(v,n=1),ytop = y_max <- par("usr")[4],col = rgb(1,0,0,0.2),border = rgb(1,0,0,0.5))
})
abline(h=-log10(0.05),col="green")
#abline(v=which(p_mapped<0.05),col=rgb(0,1,0,0.1),lwd=1)
plot(1:length(p_mapped),abs(fc[markerGroups]),type="h")
lapply(ref_hits,FUN=function(v){
  rect(xleft = v[1],ybottom = 0,xright = tail(v,n=1),ytop = y_max <- par("usr")[4],col = rgb(1,0,0,0.2),border = rgb(1,0,0,0.5))
})
plot(-log10(p_mapped),abs(fc[markerGroups]))
