#####set global parameters######################################################
setwd("/cellfile/cellnet/RF_QTL_package/")
plot_path <- "plots/swiq_exploration/"
overwrite <- F

#####load code##################################################################
source("../phosphoQTL/lib/QTL_draft.R") #has to be loaded before rfqtl
library(RFQTL)
library(ranger)
library(RColorBrewer)
source("scripts/swiq.R")


#####load data##################################################################
load("data/eQtlMappingData160728.RData")
load("data/eQTL_results160831.RData")
load("data/plotPack.RData")
names(chr_length) <- gsub(pattern = "chr",
                          replacement = "",
                          x = names(chr_length))
geno <- read.table("data/genotype_for_mapping.tsv",header=T)
chrVec <- geno[,1]
geno <- geno[,4:ncol(geno)]
geno <- t(geno)

#####results####################################################################
if(overwrite){
  #map real traits and permuted traits with min.bucket=1
  real_p_MB1 <- t(apply(phenotype[1:50,],1,
                        FUN = swiq,min.bucket=1,x=genotype[,1:3593]))
  perm_p_MB1 <- t(apply(phenotype[1:50,],1,FUN = function(y){
    idx <- sample(1:length(y))
    swiq(y=y[idx],x=genotype[,1:3593],min.bucket = 1)
  }))
  
  min(p.adjust(perm_p_MB1,method="fdr"))
  
  #map simulated traits with and without correlated predictors with min.bucket=1
  w_cor_p <- t(sapply(1:100,FUN=function(i){
    x <- matrix(NA,ncol=10,nrow=112)
    x[,1] <- sample(x = c(0,1),size = nrow(x),replace = T)
    for(i in 2:ncol(x)){
      x[,i] <- x[,i-1]
      #change only two samples for each marker
      x[sample(1:nrow(x),size = 2),i] <- sample(c(0,1),size = 2,replace = T)
    }
    y <- rnorm(n = 112)
    p <- swiq(y=y,x=x,min.bucket = 1)
    return(p)
  }))
  
  min(p.adjust(w_cor_p,method="fdr"))
  
  no_cor_p <- t(sapply(1:100,FUN=function(i){
    x <- matrix(NA,ncol=10,nrow=112)
    x[,1] <- sample(x = c(0,1),size = nrow(x),replace = T)
    for(i in 2:ncol(x)){
      x[,i] <- x[,i-1]
      #change all samples for each marker
      x[sample(1:nrow(x),size = 112),i] <- sample(c(0,1),size = 112,replace = T)
    }
    y <- rnorm(n = 112)
    p <- swiq(y=y,x=x,min.bucket = 1)
    return(p)
  }))
  
  min(p.adjust(no_cor_p,method="fdr"))
  
  #map a simulated trait with very shallow trees
  w_cor_p_shallow <- t(sapply(1:100,FUN=function(i){
    x <- matrix(NA,ncol=10,nrow=112)
    x[,1] <- sample(x = c(0,1),size = nrow(x),replace = T)
    for(i in 2:ncol(x)){
      x[,i] <- x[,i-1]
      #change only two samples for each marker
      x[sample(1:nrow(x),size = 2),i] <- sample(c(0,1),size = 2,replace = T)
    }
    y <- rnorm(n = 112)
    p <- swiq(y=y,x=x,min.bucket = 1,min.node.size = 100)
    return(p)
  }))
  
  #map a simulated trait with very shallow trees
  w_cor_p_MB15 <- t(sapply(1:100,FUN=function(i){
    x <- matrix(NA,ncol=10,nrow=112)
    x[,1] <- sample(x = c(0,1),size = nrow(x),replace = T)
    for(i in 2:ncol(x)){
      x[,i] <- x[,i-1]
      #change only two samples for each marker
      x[sample(1:nrow(x),size = 2),i] <- sample(c(0,1),size = 2,replace = T)
    }
    y <- rnorm(n = 112)
    p <- swiq(y=y,x=x,min.bucket = 15)
    return(p)
  }))
  
  #map simulated traits with imbalanced allele frequencies
  set.seed(15)
  y <- rnorm(1000)
  afreq <- rep(seq(100,500,100),each=40)
  x <- matrix(0,ncol=length(afreq),nrow=1000)
  for(i in 1:ncol(x)){
    x[sample(1:nrow(x),size = afreq[i]),i] <- 1
  }
  x <- as.data.frame(x)
  rf <- ranger(y = y,
               x = x,
               importance="permutation",
               local.importance = T,
               num.trees = 10000,
               num.threads = 0,
               min.bucket = 15)
  imps <- t(rf$variable.importance.local)
  imp_med <- t(sapply(1:ncol(x),FUN=function(i){
    m0 <- median(imps[i,x[,i]==0])
    m1 <- median(imps[i,x[,i]==1])
    return(c(m0,m1))
  }))
  imp_mean <- t(sapply(1:ncol(x),FUN=function(i){
    m0 <- mean(imps[i,x[,i]==0])
    m1 <- mean(imps[i,x[,i]==1])
    return(c(m0,m1))
  }))
  p_no_correction <- swiq(y = y,
                          x = x,
                          num.trees = 10000,
                          min.bucket = 15,
                          num.threads = 0,
                          correct_alleles = F)
  p_with_correction <- swiq(y = y,
                          x = x,
                          num.trees = 10000,
                          min.bucket = 15,
                          num.threads = 0,
                          correct_alleles = T)
  p_with_correction2 <- swiq(y = y,
                            x = x,
                            num.trees = 10000,
                            min.bucket = 15,
                            num.threads = 0,
                            correct_alleles = T,
                            correction_mode = 2,
                            window_size = 50)
  p_with_correction3 <- swiq(y = y,
                             x = x,
                             num.trees = 100000,
                             min.bucket = 15,
                             num.threads = 0,
                             correct_alleles = T,
                             correction_mode = 2)
  correction_factors <- estimate_correction(x = x,n_sim = 100)
  correction_factors2 <- estimate_correction(x = x,n_sim = 20,num.trees = 100000)
  p_with_correction4 <- swiq(y = y,
                             x = x,
                             num.trees = 10000,
                             min.bucket = 15,
                             num.threads = 0,
                             correct_alleles = T,
                             correction_mode = 3,
                             correction_mat = correction_factors)
  p_with_correction5 <- swiq(y = y,
                             x = x,
                             num.trees = 100000,
                             min.bucket = 15,
                             num.threads = 0,
                             correct_alleles = T,
                             correction_mode = 3,
                             correction_mat = correction_factors2)
  p_with_correction6 <- swiq(y = y,
                            x = x,
                            num.trees = 100000,
                            min.bucket = 15,
                            num.threads = 0,
                            correct_alleles = T)
  t10 <- swiq(y = phenotype[10,],
              x = genotype[,1:3593],
              num.trees = 10000,
              min.bucket = 15,
              correct_alleles = T,
              NA.correction = T)
  imps_cor <- correct_importances(imps = imps,x = x)
  
  set.seed(1)
  correction_factors3 <- estimate_correction(x = x,n_sim = 1)
  
  set.seed(1)
  correction_factors4 <- estimate_correction(x = x,n_sim = 1)
  
  correction_factors5 <- correction_factors
  correction_factors5[correction_factors5<0] <- 0
  correction_factors6 <- correction_factors2
  correction_factors6[correction_factors6<0] <- 0
  p_with_correction7 <- swiq(y = y,
                             x = x,
                             num.trees = 100000,
                             min.bucket = 15,
                             num.threads = 0,
                             correct_alleles = T,
                             correction_mode = 3,
                             correction_mat = correction_factors6)
  
  #compute the differences in mean between major and minor allele over
  # many traits.
  mean_deltas <- sapply(1:200,FUN=function(i){
    set.seed(i)
    y <- rnorm(1000)
    afreq <- rep(seq(100,500,100),each=40)
    x <- matrix(0,ncol=length(afreq),nrow=1000)
    for(i in 1:ncol(x)){
      x[sample(1:nrow(x),size = afreq[i]),i] <- 1
    }
    x <- as.data.frame(x)
    rf <- ranger(y = y,
                 x = x,
                 importance="permutation",
                 local.importance = T,
                 num.trees = 10000,
                 num.threads = 0,
                 min.bucket = 15)
    imps <- t(rf$variable.importance.local)
    imp_mean <- t(sapply(1:ncol(x),FUN=function(i){
      m0 <- mean(imps[i,x[,i]==0])
      m1 <- mean(imps[i,x[,i]==1])
      return(c(m0,m1))
    }))
    delta <- imp_mean[,1]-imp_mean[,2]
    return(delta)
  })
  #summarize by allele frequencies
  mean_deltas_sum <- matrix(as.vector(t(mean_deltas)),ncol=5)
  
  sim_bi_MAF <- c("p_no_correction",
                  "p_with_correction",
                  "p_with_correction2",
                  "p_with_correction3",
                  "p_with_correction4",
                  "p_with_correction5",
                  "p_with_correction6",
                  "y",
                  "x",
                  "afreq",
                  "imps",
                  "imp_med",
                  "imp_mean",
                  "t10",
                  "mean_deltas",
                  "mean_deltas_sum",
                  "imps_cor",
                  "correction_factors",
                  "correction_factors2",
                  "correction_factors3",
                  "correction_factors4")
  sim_bi_MAF <- named_list(sim_bi_MAF)
  
  #map permuted traits with real genotypes and min.bucket=15
  perm_p_MB15 <- t(apply(phenotype[1:50,],1,FUN = function(y){
    idx <- sample(1:length(y))
    swiq(y=y[idx],x=genotype[,1:3593],min.bucket = 15)
  }))
  
  #map real traits with real genotypes and min.bucket=15
  real_p_MB15 <- t(apply(phenotype[1:50,],
                         1,
                         FUN = swiq,
                         min.bucket=15,
                         x=genotype[,1:3593]))
  
  #map all traits with and without permutations
  #no perms
  i <- 0
  eQTL_p <- t(apply(phenotype,1,FUN = function(y){
    i <<- i+1
    swiq(y=y,x=genotype[,1:3593],min.bucket = 15,num.tree=10000)
  }))
  eQTL_q <- eQTL_p
  eQTL_q[] <- p.adjust(eQTL_q,method="fdr")
  
  #with perms
  i <- 0
  eQTL_p_perm <- t(apply(phenotype,1,FUN = function(y){
    i <<- i+1
    idx <- sample(1:length(y))
    swiq(y=y[idx],x=genotype[,1:3593],min.bucket = 15,num.tree=10000)
  }))
  eQTL_q_perm <- eQTL_p_perm
  eQTL_q_perm[] <- p.adjust(eQTL_q_perm,method="fdr")
  
  #with pop structure
  i <- 0
  eQTL_p_PS <- t(apply(phenotype,1,FUN = function(y){
    i <<- i+1
    swiq(y=y,x=genotype,min.bucket = 15,num.tree=10000)
  }))
  eQTL_q_PS <- eQTL_p_PS[,1:3593]
  eQTL_q_PS[] <- p.adjust(eQTL_q_PS,method="fdr")
  
  #call QTLs for the nonpermuted traits without and with modelling of pop
  #structure
  fdrs <- c(1,5,10,15,20,25) 
  eQTL <- lapply(fdrs,FUN=function(fdr){
    QTLgrouper(pmat = eQTL_q,
               sigThreshold = fdr/100,
               corThreshold = 0.8,
               distThreshold = 9,
               genotype = geno,
               chrVec = chrVec)
  })
  names(eQTL) <- paste0("FDR",fdrs)
  
  eQTL_PS <- lapply(fdrs,FUN=function(fdr){
    QTLgrouper(pmat = eQTL_q_PS,
               sigThreshold = fdr/100,
               corThreshold = 0.8,
               distThreshold = 9,
               genotype = geno,
               chrVec = chrVec)
  })
  names(eQTL_PS) <- paste0("FDR",fdrs)
  
  #not allowing NAs to inform the PI
  i <- 0
  eQTL_p_noNA <- t(apply(phenotype,1,FUN = function(y){
    i <<- i+1
    swiq(y=y,
         x=genotype[,1:3593],
         min.bucket = 15,
         num.tree=10000,
         NA.correction = T,
         multi=T)
  }))
  eQTL_q_noNA <- eQTL_p_noNA
  eQTL_q_noNA[] <- p.adjust(eQTL_q_noNA,method="fdr")
  
  i <- 0
  eQTL_p_perm_noNA <- t(apply(phenotype,1,FUN = function(y){
    i <<- i+1
    idx <- sample(1:length(y))
    swiq(y=y[idx],
         x=genotype[,1:3593],
         min.bucket = 15,
         num.tree=10000,
         NA.correction = T,
         multi=T)
  }))
  eQTL_q_perm_noNA <- eQTL_p_perm_noNA
  eQTL_q_perm_noNA[] <- p.adjust(eQTL_q_perm_noNA,method="fdr")
  
  eQTL_noNA <- lapply(fdrs,FUN=function(fdr){
    QTLgrouper(pmat = eQTL_q_noNA,
               sigThreshold = fdr/100,
               corThreshold = 0.8,
               distThreshold = 9,
               genotype = geno,
               chrVec = chrVec)
  })
  names(eQTL_noNA) <- paste0("FDR",fdrs)
  
  #use a t-test instead of the binomial
  i <- 0
  eQTL_p_noNA_t <- t(apply(phenotype,1,FUN = function(y){
    i <<- i+1
    swiq(y=y,
         x=genotype[,1:3593],
         min.bucket = 15,
         num.tree=10000,
         NA.correction = T,
         multi=T,
         method = "t-test")
  }))
  eQTL_q_noNA_t <- eQTL_p_noNA_t
  eQTL_q_noNA_t[] <- p.adjust(eQTL_q_noNA_t,method="fdr")
  eQTL_noNA_t <- lapply(fdrs,FUN=function(fdr){
    QTLgrouper(pmat = eQTL_q_noNA_t,
               sigThreshold = fdr/100,
               corThreshold = 0.8,
               distThreshold = 9,
               genotype = geno,
               chrVec = chrVec)
  })
  names(eQTL_noNA_t) <- paste0("FDR",fdrs)
  
  i <- 0
  eQTL_p_perm_noNA_t <- t(apply(phenotype,1,FUN = function(y){
    i <<- i+1
    idx <- sample(1:length(y))
    swiq(y=y[idx],
         x=genotype[,1:3593],
         min.bucket = 15,
         num.tree=10000,
         NA.correction = T,
         multi=T,
         method = "t-test")
  }))
  eQTL_q_perm_noNA_t <- eQTL_p_perm_noNA_t
  eQTL_q_perm_noNA_t[] <- p.adjust(eQTL_q_perm_noNA_t,method="fdr")
  
  save(real_p_MB1,
       perm_p_MB1,
       real_p_MB15,
       perm_p_MB15,
       w_cor_p,
       no_cor_p,
       w_cor_p_shallow,
       eQTL_p,
       eQTL_q,
       eQTL_p_perm,
       eQTL_q_perm,
       eQTL_p_PS,
       eQTL_q_PS,
       eQTL,
       eQTL_PS,
       eQTL_p_noNA,
       eQTL_q_noNA,
       eQTL_noNA,
       eQTL_p_perm_noNA,
       eQTL_q_perm_noNA,
       eQTL_p_perm_noNA_t,
       eQTL_q_perm_noNA_t,
       eQTL_p_noNA_t,
       eQTL_q_noNA_t,
       eQTL_noNA_t,
       sim_bi_MAF,
       file="data/swiq_exploration.RData")
}else{
  load("data/swiq_exploration.RData")
}


#####plots######################################################################
# qtl plots for traits 7 and 20 with the old and new (no min bucket) approaches
# trait 10
png(paste0(plot_path,"profiles_MB1_trait10.png"),width = 700,height=700)
par(mfrow=c(2,1))
QTL_profile(p_vec = real_p_MB1[10,],
            predictorLocs = markerLocs,
            chrLen = chr_length,
            xlab="marker position",
            chrCol = "red",
            ylab = "-log10(p)",
            main = "SWIQ approach",
            col = "black",
            chrLty = 2,
            labcex = 1.5)

QTL_profile(p_vec = pv[10,],
            predictorLocs = markerLocs,
            chrLen = chr_length,
            xlab="marker position",
            chrCol = "red",
            ylab = "-log10(p)",
            main = "permutation approach",
            col = "black",
            chrLty = 2,
            labcex = 1.5)
dev.off()

# trait 26
png(paste0(plot_path,"profiles_MB1_trait26.png"),width = 700,height=700)
par(mfrow=c(2,1))
QTL_profile(p_vec = real_p_MB1[26,],
            predictorLocs = markerLocs,
            chrLen = chr_length,
            xlab="marker position",
            chrCol = "red",
            ylab = "-log10(p)",
            main = "SWIQ approach",
            col = "black",
            chrLty = 2,
            labcex = 1.5)

QTL_profile(p_vec = pv[26,],
            predictorLocs = markerLocs,
            chrLen = chr_length,
            xlab="marker position",
            chrCol = "red",
            ylab = "-log10(p)",
            main = "permutation approach",
            col = "black",
            chrLty = 2,
            labcex = 1.5)
dev.off()

#histogramm of p-values for non correlated predictors
p <- gg_histogram(x = as.vector(no_cor_p),
                  breaks=10,
                  title = "p-values for uncorrelated predictors",
                  xlab = "p-value")
ggsave(filename = paste0(plot_path,"no_cor_p_hist.png"),scale = 2)

#histogramm of p-values for correlated predictors
p <- gg_histogram(x = as.vector(w_cor_p),
                  breaks=10,
                  title = "p-values for correlated predictors",
                  xlab = "p-value")
ggsave(filename = paste0(plot_path,"w_cor_p_hist.png"),scale = 2)

#histogramm of p-values for correlated predictors
p <- gg_histogram(x = as.vector(w_cor_p_shallow),
                  breaks=10,
                  title = "p-values for correlated predictors, shallow trees",
                  xlab = "p-value")
ggsave(filename = paste0(plot_path,"w_cor_shallow_p_hist.png"),scale = 2)

p <- gg_histogram(x = as.vector(w_cor_p_MB15),
                  breaks=10,
                  title = "p-values for correlated predictors, no small nodes",
                  xlab = "p-value")
ggsave(filename = paste0(plot_path,"w_cor_MB15_p_hist.png"),scale = 2)

#QTL map at FDR<10% with old results
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]
names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")
png(filename = paste0(plot_path,"qtlmap_legacy.png"),width = 960,height = 960)
par(cex=1.5,cex.lab=2,mar=c(5,5,4,2)+0.1)
QTLplotter(QTLlist = QtlList$FDR10,
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
           main="permutation-based results")
dev.off()

#new results
png(filename = paste0(plot_path,"qtlmap_swiq.png"),width = 960,height = 960)
par(cex=1.5,cex.lab=2,mar=c(5,5,4,2)+0.1)
QTLplotter(QTLlist = eQTL$FDR10,
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
           main="swiq results")
dev.off()

#new results
png(filename = paste0(plot_path,"qtlmap_swiq_fdr1.png"),
    width = 960,
    height = 960)
par(cex=1.5,cex.lab=2,mar=c(5,5,4,2)+0.1)
QTLplotter(QTLlist = eQTL$FDR1,
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
           main="swiq results")
dev.off()

#plot only lead markers
eQTL_lead <- lapply(eQTL$FDR10,FUN=function(qtl){
  qtl$predictors <- matrix(rep(qtl$mostSignificantPredictor,2),nrow=1,ncol=2)
  colnames(qtl$predictors) <- c("start","end")
  return(qtl)
})

png(filename = paste0(plot_path,"qtlmap_swiq_lead_only.png"),
    width = 960,
    height = 960)
par(cex=1.5,cex.lab=2,mar=c(5,5,4,2)+0.1)
QTLplotter(QTLlist = eQTL_lead,
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
           main="swiq results, FDR<10%, only lead marker")
dev.off()

eQTL_lead_fdr1 <- lapply(eQTL$FDR1,FUN=function(qtl){
  qtl$predictors <- matrix(rep(qtl$mostSignificantPredictor,2),nrow=1,ncol=2)
  colnames(qtl$predictors) <- c("start","end")
  return(qtl)
})

png(filename = paste0(plot_path,"qtlmap_swiq_lead_only_fdr1.png"),
    width = 960,
    height = 960)
par(cex=1.5,cex.lab=2,mar=c(5,5,4,2)+0.1)
QTLplotter(QTLlist = eQTL_lead_fdr1,
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
           main="swiq results, FDR<1%, only lead marker")
dev.off()

#is there a relationship between amount of NAs and p-values?
gg_histogram(colSums(eQTL_p_perm_noNA<0.1),title = "")
gg_histogram(colSums(eQTL_p_perm<0.1),title = "")
gg_barplot(colSums(eQTL_p_perm<0.1),title = "")
p <- gg_boxplot(x=colSums(is.na(genotype[,1:3593]))>0,
                y=colSums(eQTL_p_perm<0.1),
                xlab = "NAs in the predictor",
                ylab = "sum of p-values under 0.1",
                title="no correction for NAs")
ggsave(filename = paste0(plot_path,"boxplot_perms_no_NA_correction.png"),
       scale = 2)
p <- gg_boxplot(x=colSums(is.na(genotype[,1:3593]))>0,
                y=colSums(eQTL_p_perm_noNA<0.1),
                xlab = "NAs in the predictor",
                ylab = "sum of p-values under 0.1",
                title="with correction for NAs")
ggsave(filename = paste0(plot_path,"boxplot_perms_with_NA_correction.png"),
       scale = 2)
p <- gg_boxplot(x=colSums(is.na(genotype[,1:3593]))>0,
                y=colSums(eQTL_p_perm_noNA<0.01),
                xlab = "NAs in the predictor",
                ylab = "sum of p-values under 0.01",
                title="with correction for NAs")
ggsave(filename = paste0(plot_path,"boxplot_perms_with_NA_correction_001.png"),
       scale = 2)

eQTL_lead_fdr10_NA_correction <- lapply(eQTL_noNA$FDR10,FUN=function(qtl){
  qtl$predictors <- matrix(rep(qtl$mostSignificantPredictor,2),nrow=1,ncol=2)
  colnames(qtl$predictors) <- c("start","end")
  return(qtl)
})

png(filename = paste0(plot_path,
                      "qtlmap_swiq_lead_only_fdr10_NA_correction.png"),
    width = 960,
    height = 960)
par(cex=1.5,cex.lab=2,mar=c(5,5,4,2)+0.1)
QTLplotter(QTLlist = eQTL_lead_fdr10_NA_correction,
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
           main="swiq results, FDR<10%, only lead marker, NA corrected")
dev.off()

eQTL_lead_fdr10_NA_correction_t <- lapply(eQTL_noNA_t$FDR10,FUN=function(qtl){
  qtl$predictors <- matrix(rep(qtl$mostSignificantPredictor,2),nrow=1,ncol=2)
  colnames(qtl$predictors) <- c("start","end")
  return(qtl)
})

png(filename = paste0(plot_path,
                      "qtlmap_swiq_lead_only_fdr10_NA_correction_t.png"),
    width = 960,
    height = 960)
par(cex=1.5,cex.lab=2,mar=c(5,5,4,2)+0.1)
QTLplotter(QTLlist = eQTL_lead_fdr10_NA_correction_t,
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
           main="swiq results, FDR<10%, t-test")
dev.off()

#make plots for approach to correct for imbalanced allele frequencies
png(filename = paste0(plot_path,
                      "med_comp.png"),
    width = 600,
    height = 600)
plot(sim_bi_MAF$imp_med[,1:2],
     ylab="median of minor allele",
     xlab="median of major allele",
     pch=19,
     col=rgb(1-sim_bi_MAF$afreq/500,0,0,0.7))
abline(h=0,v=0)
dev.off()
cor.test(sim_bi_MAF$imp_med[,1],sim_bi_MAF$imp_med[,2],method="spearman")
gg_histogram(sim_bi_MAF$p_no_correction,xlab="p-value",
             title="no correction for allele imbalances")
ggsave(filename = paste0(plot_path,"histogram_p_MAF_bias.png"),
       scale = 2)
gg_histogram(sim_bi_MAF$p_with_correction,xlab="p-value",
             title="with correction for allele imbalances")
ggsave(filename = paste0(plot_path,"histogram_p_no_MAF_bias.png"),
       scale = 2)
gg_boxplot(y=sim_bi_MAF$p_no_correction,x=sim_bi_MAF$afreq,
           xlab = "#of samples with minor allele out of 1000",ylab="p-value")
ggsave(filename = paste0(plot_path,"boxplot_p_MAF_bias.png"),
       scale = 2)
m <- which.min(sim_bi_MAF$p_no_correction)
gg_boxplot(y=sim_bi_MAF$imps[m,],
           x=sim_bi_MAF$x[,m],
           box_labels = paste0(c(0,1),", n = ",c(sum(sim_bi_MAF$x[,m]==0),sum(sim_bi_MAF$x[,m]==1))),
           xlab="alleles",
           ylab="sample wise permutation importances",
           title="raw importance values")
ggsave(filename = paste0(plot_path,"boxplot_imps_example.png"),
       scale = 2)
gg_boxplot(y=sim_bi_MAF$imps_cor[m,],
           x=sim_bi_MAF$x[,m],
           box_labels = paste0(c(0,1),", n = ",c(sum(sim_bi_MAF$x[,m]==0),sum(sim_bi_MAF$x[,m]==1))),
           xlab="alleles",
           ylab="sample wise permutation importances",
           title="corrected importance values")
ggsave(filename = paste0(plot_path,"boxplot_imps_corrected_example.png"),
       scale = 2)
m <- which(imp_med[,1]>0&imp_med[,2]>0)[1]
gg_boxplot(y=sim_bi_MAF$imps[m,],
           x=sim_bi_MAF$x[,m],
           box_labels = paste0(c(0,1),", n = ",c(sum(sim_bi_MAF$x[,m]==0),sum(sim_bi_MAF$x[,m]==1))),
           xlab="alleles",
           ylab="sample wise permutation importances",
           title="raw importance values")
ggsave(filename = paste0(plot_path,"boxplot_imps_example2.png"),
       scale = 2)
gg_boxplot(y=sim_bi_MAF$imps_cor[m,],
           x=sim_bi_MAF$x[,m],
           box_labels = paste0(c(0,1),", n = ",c(sum(sim_bi_MAF$x[,m]==0),sum(sim_bi_MAF$x[,m]==1))),
           xlab="alleles",
           ylab="sample wise permutation importances",
           title="corrected importance values")
ggsave(filename = paste0(plot_path,"boxplot_imps_corrected_example2.png"),
       scale = 2)
png(filename = paste0(plot_path,
                      "trait_10_corrected.png"),
    width = 600,
    height = 600)
plot(-log10(sim_bi_MAF$t10),-log10(eQTL_p_noNA[10,]),xlab="-log10(p) after correction",ylab="-log10(p) without correction")
abline(a=0,b=1)
dev.off()

# show the difference of the means by allele and MAF
gg_boxplot(y=sim_bi_MAF$imp_mean[,1]-sim_bi_MAF$imp_mean[,2],
           x=sim_bi_MAF$afreq/10,
           xlab="MAF",
           ylab="(mean SWIPs major allele)-(mean SWIPs minor allele)",
           title="SWPI difference between alleles")
ggsave(filename = paste0(plot_path,"boxplot_delta_mean.png"),
       scale = 2)

gg_boxplot(y=as.vector(sim_bi_MAF$mean_deltas_sum),
           x=rep(c(10,20,30,40,50),each=8000),
           xlab="MAF",
           ylab="(mean SWIPs major allele)-(mean SWIPs minor allele)",
           title="SWPI difference between alleles, 200 traits, no outliers",
           show_outliers = F)
ggsave(filename = paste0(plot_path,"boxplot_delta_sum_mean.png"),
       scale = 2)

#show pvalues after correction with sliding window
gg_boxplot(y=sim_bi_MAF$p_with_correction2,x=sim_bi_MAF$afreq,
           xlab = "#of samples with minor allele out of 1000",ylab="p-value")
ggsave(filename = paste0(plot_path,"boxplot_p_sliding_window.png"),
       scale = 2)

#show p-values after correction with 10k and 100k trees
gg_boxplot(y=sim_bi_MAF$p_with_correction4,x=sim_bi_MAF$afreq,
           xlab = "#of samples with minor allele out of 1000",
           ylab="p-value",
           title = "estimated correction factors, 10k trees")
ggsave(filename = paste0(plot_path,"boxplot_p_est_10k.png"),
       scale = 2)

gg_boxplot(y=sim_bi_MAF$p_with_correction5,x=sim_bi_MAF$afreq,
           xlab = "#of samples with minor allele out of 1000",
           ylab="p-value",
           title = "estimated correction factors, 100k trees")
ggsave(filename = paste0(plot_path,"boxplot_p_est_100k.png"),
       scale = 2)

#show that the allele bias is trait dependent 
gg_dotplot(x = correction_factors3[afreq==100,1],
           y = correction_factors4[afreq==100,1],
           pch=19,
           main="major allele SWPIs, same trait, markers with MAF=10%",
           xlab = "simulation 1",
           ylab = "simulation 2")
ggsave(filename = paste0(plot_path,"dotplot_est_twice.png"),
       scale = 2)
