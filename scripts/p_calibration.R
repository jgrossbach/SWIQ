#####set global parameters######################################################
setwd("/cellfile/cellnet/RF_QTL_package/")
plot_path <- "plots/swiq_exploration/"
eQTL_path <- "data/mapped_albert_traits_small_sample/"
overwrite <- F

#####load code##################################################################
source("../phosphoQTL/lib/QTL_draft.R") #has to be loaded before rfqtl
library(RFQTL)
library(ranger)
library(RColorBrewer)
library(biomaRt)
source("scripts/swiq.R")
source("scripts/help_fun.R")


#####load data##################################################################

albert_eqtl <- readRDS("data/eQTL_albert_results.rds")
mapping_data <- readRDS(file = "data/albert_mapping_data.rds")
genotype <- mapping_data$genotype
phenotype <- mapping_data$phenotype
rownames(phenotype) <- gsub("[.]","-",rownames(phenotype))
folds <- readRDS(file = "data/albert_folds.rds")
chrVec <- colnames(mapping_data$genotype_full)
chrVec <- gsub(pattern = "[.].*","",x = chrVec)
chrVec <- chrVec[sapply(1:ncol(genotype),FUN=function(m){which(mapping_data$marker_groups==m)[1]})]
snp_loc <- sapply(colnames(mapping_data$genotype_full),FUN=function(x){
  pos <- gsub(pattern = "^[^.]*[.]","",x)
  pos <- gsub(pattern = "_.*","",pos)
  return(as.numeric(pos))
})
m_loc <- t(sapply(1:ncol(genotype),FUN=function(m){
  range(snp_loc[mapping_data$marker_groups==m])
}))
colnames(m_loc) <- c("start","end")
m_loc <- data.frame(chr=chrVec,m_loc,m_mid=ceiling(rowMeans(m_loc)))

load("data/plotPack.RData")
names(chr_length) <- gsub(pattern = "chr",
                          replacement = "",
                          x = names(chr_length))

geno <- read.table("data/genotype_for_mapping.tsv",header=T)
#chrVec <- geno[,1]
geno <- geno[,4:ncol(geno)]
geno <- t(geno)
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]
names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")

# geno <- read.table("data/genotype_for_mapping.tsv",header=T)
# chrVec <- geno[,1]
# mloc <- geno[,1:3]
# geno <- geno[,4:ncol(geno)]
# geno <- t(geno)

#biomart
mart <- useMart("ensembl")
ensembl <- useDataset("scerevisiae_gene_ensembl", mart = mart)


eQTL_phospho <- read.table("data/eQTL.tsv",sep="\t",header=T)
#####results####################################################################


loc_fcs <- t(sapply(1:nrow(phenotype),FUN=function(trait){
  trait <<- trait
  set.seed(trait)
  s <- sample(1:ncol(phenotype),size = 200,replace = F)
  phenotype_s <- phenotype[trait,s]
  genotype_s <- genotype[s,]
  fc <- sapply(1:ncol(genotype),FUN=function(m){
    mean(phenotype_s[genotype_s[,m]==1])-mean(phenotype_s[genotype_s[,m]==0])
  })
  return(fc)
}))
loc_fcs2 <- t(sapply(1:nrow(phenotype),FUN=function(trait){
  trait <<- trait
  set.seed(trait)
  s <- sample(1:ncol(phenotype),size = 200,replace = F)
  phenotype_s <- phenotype[trait,-s]
  genotype_s <- genotype[-s,]
  fc <- sapply(1:ncol(genotype),FUN=function(m){
    mean(phenotype_s[genotype_s[,m]==1])-mean(phenotype_s[genotype_s[,m]==0])
  })
  return(fc)
}))

same_direction_m <- (loc_fcs>0&loc_fcs2>0)|(!loc_fcs>0&!loc_fcs2>0)

full_fc <- sapply(1:ncol(genotype),FUN=function(m){
  m <<- m
  rowMeans(phenotype[,genotype[,m]==1],na.rm=T)-rowMeans(phenotype[,genotype[,m]==0],na.rm=T)
})

same_direction_m2 <- (full_fc>0&loc_fcs2>0)|(!full_fc>0&!loc_fcs2>0)
same_direction_m3 <- (full_fc>0&loc_fcs>0)|(!full_fc>0&!loc_fcs>0)

adjusted_fc <- (loc_fcs*200+full_fc*(ncol(phenotype)-200))/ncol(phenotype)

files <- list.files(eQTL_path)
present_files <- paste0("trait_",1:nrow(phenotype),".rds")%in%files
last_file <- which(!present_files)[1]-1
if(is.na(last_file)){
  last_file <- length(present_files)
}
p <- t(sapply(1:last_file,FUN=function(i){
  fname <- paste0(eQTL_path,"trait_",i,".rds")
  out <- readRDS(fname)
  return(out[,1])
}))

max_out <- apply(phenotype,1,FUN=function(x){
  max(abs(x-mean(x)))
})
#p_clean <- which(max_out>quantile(max_out,probs = 0.5))
#p_clean <- which(n_out>20)
p_clean <- 1:nrow(p)

q <- p[p_clean,]
q[] <- p.adjust(q,method="fdr")

qvec <- as.vector(q)
fc <- as.vector(loc_fcs[p_clean,])
# fc2 <- as.vector(loc_fcs2[1:last_file,])
#fc2 <- as.vector(folds$fcs[[setdiff(1:2,i)]][1:last_file,])
# fc2 <- (fc+fc2)/2
# fc <- random_fc[1:nrow(qmat),]
# fc2 <- random_fc2[1:nrow(qmat),]
#fc2 <- as.vector(loc_fcs2[p_clean,])
fc2 <- as.vector(full_fc[p_clean,])
fc <- fc>0
fc2 <- fc2>0
same_direction <- as.vector((fc&fc2)|(!fc&!fc2))
q_order <- order(qvec,sample(1:length(qvec),size=length(qvec)))
diff_cum <- cumsum(!same_direction[q_order])
fp_cum <- cumsum(qvec[q_order])
q_order_traits <- q_order%%length(p_clean)
q_order_traits[q_order_traits==0] <- length(p_clean)
fdr_mat <- data.frame(exp=fp_cum,obs=diff_cum,trait=q_order_traits,q=qvec[q_order],same_direction=same_direction[q_order])
fdr_sub <- fdr_mat[which(!fdr_mat[,"same_direction"]&fdr_mat[,4]<0.01),]
fdr_sub2 <- fdr_sub[fdr_sub$trait<=nrow(albert_eqtl$q),]
#wrong_dir <- which(q<0.01&!same_direction_m,arr.ind=T)

#par(mfrow=c(1,1))
plot(fdr_mat[1:sum(q<0.2),1:2])
abline(a=0,b=1,col="red")
abline(a=0,b=0.5,col="green")

plot()

plot(fdr_mat[1:sum(q<0.01),])
abline(a=0,b=1,col="red")
abline(a=0,b=0.5,col="green")

fp_m <- colSums(q<0.01&!same_direction_m)
fp_t <- rowSums(q<0.01&!same_direction_m)

m <- 2543
q_thr <- 0.05

fp_slopes <- sapply(1:ncol(genotype),FUN=function(m){
  m <<- m
  p_m <- p[,m]
  q_m <- p.adjust(p = p_m,method = "fdr")
  o_m <- order(p_m,sample(1:nrow(p)))
  if(sum(q_m<=q_thr)<10){return(NA)}
  o_m <- o_m[1:sum(q_m<=q_thr)]

  fp_lm <- lm(cumsum(!same_direction_m[o_m,m])~cumsum(q_m[o_m]))
  return(fp_lm$coefficients[2])
})

fp_slopes2 <- sapply(1:ncol(genotype),FUN=function(m){
  m <<- m
  p_m <- p[,m]
  q_m <- p.adjust(p = p_m,method = "fdr")
  o_m <- order(p_m,sample(1:nrow(p)))
  #if(sum(q_m<=q_thr)<10){return(NA)}
  o_m <- o_m[1:1000]
  
  fp_lm <- lm(cumsum(!same_direction_m[o_m,m])~cumsum(q_m[o_m]))
  return(fp_lm$coefficients[2])
})

m <- 2545
p_m <- p[,m]
q_m <- p.adjust(p = p_m,method = "fdr")
o_m <- order(p_m,sample(1:nrow(p)))
#o_m <- o_m[1:sum(q_m<=q_thr)]
o_m <- o_m[1:1000]
plot(cumsum(q_m[o_m]),cumsum(!same_direction_m[o_m,m]))
abline(a=0,b=1,col="red")
abline(a=0,b=0.5,col="green")

fdr_plots2 <- lapply(1:2,FUN=function(i){
  qmat <- q[[i]]
  qfull <- albert_eqtl$q
  nt <- min(nrow(qmat),nrow(qfull))
  qmat <- qmat[1:nt,]
  qfull <- qfull[1:nt,]
  qvec <- as.vector(qmat)
  q_order <- order(qvec,sample(1:length(qvec),size=length(qvec)))
  diff_cum <- cumsum(qfull[q_order]>0.1)
  return(data.frame(qvec[q_order],diff_cum))
})

fdr10_check <- t(sapply(1:last_file,FUN=function(i){
  q <- p.adjust(p[[1]][i,],method="fdr")
  hits <- q<0.1
  fc <- folds$fcs[[1]][i,hits]
  fc2 <- folds$fcs[[2]][i,hits]
  # fc2 <- (fc+fc2)/2
  # fc <- random_fc[i,hits]
  # fc2 <- random_fc2[i,hits]
  wrong_dir <- sum((fc*fc2)<0)
  return(c(sum(hits)*0.1,wrong_dir))
}))

loc_fcs <- t(sapply(1:nrow(phenotype),FUN=function(trait){
  set.seed(trait)
  s <- sample(1:ncol(phenotype),size = 200,replace = F)
  phenotype_s <- phenotype[trait,s]
  genotype_s <- genotype[s,]
  fc <- sapply(1:ncol(genotype),FUN=function(m){
    mean(phenotype_s[genotype_s[,m]==1])-mean(phenotype_s[genotype_s[,m]==0])
  })
  return(fc)
}))
