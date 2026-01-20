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
source("scripts/help_fun.R")


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
  sim_list <- list()
  
  ###small simulations###
  n_markers <- 100
  n_samples <- 200
  n_sims <- 20
  # sim1: do i find qtl for traits that should have one?
  n_hits <- 2
  effect_size <- 1
  sim_cur <- lapply(1:n_sims,FUN=function(i){
    set.seed(i)
    i <<- i
    y <- rnorm(n_samples)
    afreq <- rep(0.5,n_markers)
    x <- matrix(0,ncol=length(afreq),nrow=n_samples)
    for(j in 1:ncol(x)){
      x[sample(1:nrow(x),size = floor(afreq[j]*n_samples)),j] <- 1
    }
    hits <- sample(1:n_markers,size = n_hits,replace = F)
    for(j in 1:n_hits){
      m <- hits[j]
      affected <- x[,m]==1
      if(j%%2==0){
        y[affected] <- y[affected]+effect_size
      }else{
        y[affected] <- y[affected]-effect_size
      }
    }
    hits_out <- rep(0,n_markers)
    for(j in 1:n_hits){
      hits_out[hits[j]] <- j
    }
    idx <- sample(1:n_samples,size = n_samples)
    p1 <- swiq(y = y,
               x = x,
               num.trees = 5000,
               num.threads = 20,
               min.bucket = 1,
               imp.blocks = 1)
    p2 <- swiq(y = y[idx],
               x = x,
               num.trees = 5000,
               num.threads = 20,
               min.bucket = 1,
               imp.blocks = 1)
    return(cbind(hits_out,p1,p2))
  })
  sim_list$sim1 <- sim_cur
  
  # sim2: add correlations to sim1
  n_hits <- 2
  effect_size <- 1
  n_replace <- floor(n_samples*0.1)
  sim_cur <- lapply(1:n_sims,FUN=function(i){
    set.seed(i)
    i <<- i
    y <- rnorm(n_samples)
    x <- matrix(0,ncol=n_markers,nrow=n_samples)
    x[,1] <- sample(x = c(0,1),size = nrow(x),replace = T)
    for(i in 2:ncol(x)){
      x[,i] <- x[,i-1]
      #change only some samples for each marker
      x[sample(1:nrow(x),size = n_replace),i] <- sample(c(0,1),
                                                        size = n_replace,
                                                        replace = T)
    }
    hits <- sample(1:n_markers,size = n_hits,replace = F)
    for(j in 1:n_hits){
      m <- hits[j]
      affected <- x[,m]==1
      if(j%%2==0){
        y[affected] <- y[affected]+effect_size
      }else{
        y[affected] <- y[affected]-effect_size
      }
    }
    hits_out <- rep(0,n_markers)
    for(j in 1:n_hits){
      hits_out[hits[j]] <- j
    }
    idx <- sample(1:n_samples,size = n_samples)
    p1 <- swiq(y = y,
               x = x,
               num.trees = 5000,
               num.threads = 20,
               min.bucket = 1,
               imp.blocks = 1)
    p2 <- swiq(y = y[idx],
               x = x,
               num.trees = 5000,
               num.threads = 20,
               min.bucket = 1,
               imp.blocks = 1)
    return(cbind(hits_out,p1,p2))
  })
  sim_list$sim2 <- sim_cur
  
  # sim3: add different allele frequencies to sim 1
  n_hits <- 2
  effect_size <- 1
  sim_cur <- lapply(1:n_sims,FUN=function(i){
    set.seed(i)
    i <<- i
    y <- rnorm(n_samples)
    afreq <- rep(x=seq(0.1,0.5,0.1),length.out = n_markers)
    afreq <- sort(afreq)
    x <- matrix(0,ncol=length(afreq),nrow=n_samples)
    for(j in 1:ncol(x)){
      x[sample(1:nrow(x),size = floor(afreq[j]*n_samples)),j] <- 1
    }
    hits <- sample(1:n_markers,size = n_hits,replace = F)
    for(j in 1:n_hits){
      m <- hits[j]
      affected <- x[,m]==1
      if(j%%2==0){
        y[affected] <- y[affected]+effect_size
      }else{
        y[affected] <- y[affected]-effect_size
      }
    }
    hits_out <- rep(0,n_markers)
    for(j in 1:n_hits){
      hits_out[hits[j]] <- j
    }
    idx <- sample(1:n_samples,size = n_samples)
    p1 <- swiq(y = y,
               x = x,
               num.trees = 5000,
               num.threads = 20,
               min.bucket = 1,
               imp.blocks = 1)
    p2 <- swiq(y = y[idx],
               x = x,
               num.trees = 5000,
               num.threads = 20,
               min.bucket = 1,
               imp.blocks = 1)
    return(cbind(hits_out,p1,p2,afreq))
  })
  sim_list$sim3 <- sim_cur
  
  # sim4: change min bucket from sim1
  n_hits <- 2
  effect_size <- 1
  sim_cur <- lapply(1:n_sims,FUN=function(i){
    set.seed(i)
    i <<- i
    y <- rnorm(n_samples)
    afreq <- rep(0.5,n_markers)
    x <- matrix(0,ncol=length(afreq),nrow=n_samples)
    for(j in 1:ncol(x)){
      x[sample(1:nrow(x),size = floor(afreq[j]*n_samples)),j] <- 1
    }
    hits <- sample(1:n_markers,size = n_hits,replace = F)
    for(j in 1:n_hits){
      m <- hits[j]
      affected <- x[,m]==1
      if(j%%2==0){
        y[affected] <- y[affected]+effect_size
      }else{
        y[affected] <- y[affected]-effect_size
      }
    }
    hits_out <- rep(0,n_markers)
    for(j in 1:n_hits){
      hits_out[hits[j]] <- j
    }
    idx <- sample(1:n_samples,size = n_samples)
    p1 <- swiq(y = y,
               x = x,
               num.trees = 5000,
               num.threads = 20,
               min.bucket = 15,
               imp.blocks = 1)
    p2 <- swiq(y = y[idx],
               x = x,
               num.trees = 5000,
               num.threads = 20,
               min.bucket = 15,
               imp.blocks = 1)
    return(cbind(hits_out,p1,p2))
  })
  sim_list$sim4 <- sim_cur
  
  #sim5: have more than two alleles
  n_hits <- 2
  effect_size <- 1
  sim_cur <- lapply(1:n_sims,FUN=function(i){
    set.seed(i)
    i <<- i
    y <- rnorm(n_samples)
    x <- matrix(0,ncol=n_markers,nrow=n_samples)
    for(j in 1:ncol(x)){
      x[,j] <- sample(c(0,1,2),size = n_samples,replace=T)
    }
    hits <- sample(1:n_markers,size = n_hits,replace = F)
    for(j in 1:n_hits){
      m <- hits[j]
      affected <- x[,m]==1
      if(j%%2==0){
        y[affected] <- y[affected]+effect_size
      }else{
        y[affected] <- y[affected]-effect_size
      }
    }
    hits_out <- rep(0,n_markers)
    for(j in 1:n_hits){
      hits_out[hits[j]] <- j
    }
    idx <- sample(1:n_samples,size = n_samples)
    p1 <- swiq(y = y,
               x = x,
               num.trees = 5000,
               num.threads = 20,
               min.bucket = 1,
               imp.blocks = 1)
    p2 <- swiq(y = y[idx],
               x = x,
               num.trees = 5000,
               num.threads = 20,
               min.bucket = 1,
               imp.blocks = 1)
    return(cbind(hits_out,p1,p2))
  })
  sim_list$sim5 <- sim_cur
  
  ###map real traits###
  #real traits were mapped with an external script
  system("Rscript scripts/swiq_eQTL.R")
  
  save(sim_list,
       file="data/swiq_exploration.RData")
}else{
  load("data/swiq_exploration.RData")
}


#####plots######################################################################
#sim1
p_list <- lapply(1:length(sim_list$sim1),FUN=function(i){
  plot_mat <- sim_list$sim1[[i]]
  max_p <- max(-log10(plot_mat[,2:3]))
  hits <- which(plot_mat[,1]!=0)
  colvec <- ifelse(1:nrow(plot_mat)%in%hits,"red","black")
  gg_dotplot(x = -log10(plot_mat[,2]),
             y = -log10(plot_mat[,3]),
             ylim=c(0,max_p),
             xlim=c(0,max_p),
             col=colvec,
             xlab = "-log10(p), non permuted",
             ylab = "-log10(p), permuted",
             main = paste0("simulation ",i))
  
})

pvec_unperm <- unlist(lapply(sim_list$sim1,FUN=function(plot_mat){
  return(plot_mat[plot_mat[,1]==0,2])
}))
gg_histogram(x = pvec_unperm,
             breaks = 20,
             color = "darkgreen",
             xlab = "p-value",
             ylab = "frequency",
             title = "non permuted traits, markers without simulated effects")
pvec_perm <- unlist(lapply(sim_list$sim1,FUN=function(plot_mat){
  return(plot_mat[plot_mat[,1]==0,3])
}))

gg_histogram(x = pvec_perm,
             breaks = 20,
             color = "darkgreen",
             xlab = "p-value",
             ylab = "frequency",
             title = "permuted traits, markers without simulated effects")

gg_boxplot(y = c(pvec_perm,pvec_unperm),
           x = rep(x = c("permuted","unpermuted"),each=length(pvec_perm)))
wilcox.test(pvec_perm,pvec_unperm)
#sim2
p_list <- lapply(1:length(sim_list$sim2),FUN=function(i){
  plot_mat <- sim_list$sim2[[i]]
  max_p <- max(-log10(plot_mat[,2:3]))
  hits <- which(plot_mat[,1]!=0)
  colvec <- ifelse(1:nrow(plot_mat)%in%hits,"red","black")
  gg_dotplot(x = -log10(plot_mat[,2]),
             y = -log10(plot_mat[,3]),
             ylim=c(0,max_p),
             xlim=c(0,max_p),
             col=colvec,
             xlab = "-log10(p), non permuted",
             ylab = "-log10(p), permuted",
             main = paste0("simulation ",i))
  
})

pvec_perm <- unlist(lapply(sim_list$sim2,FUN=function(plot_mat){
  return(plot_mat[plot_mat[,1]==0,3])
}))

gg_histogram(x = pvec_perm,
             breaks = 20,
             color = "darkgreen",
             xlab = "p-value",
             ylab = "frequency",
             title = "permuted traits, correlated markers without simulated effects")

#sim3
p_list <- lapply(1:length(sim_list$sim3),FUN=function(i){
  plot_mat <- sim_list$sim3[[i]]
  max_p <- max(-log10(plot_mat[,2:3]))
  hits <- which(plot_mat[,1]!=0)
  gg_boxplot(x = plot_mat[,4],
             y = plot_mat[,3],
             title = paste0("simulation ",i),
             xlab = "MAF",
             ylab = "p-value")
})
p_list <- lapply(1:length(sim_list$sim3),FUN=function(i){
  plot_mat <- sim_list$sim3[[i]]
  max_p <- max(-log10(plot_mat[,2:3]))
  hits <- which(plot_mat[,1]!=0)
  gg_boxplot(x = plot_mat[,1],
             y = plot_mat[,2],
             title = paste0("simulation ",i),
             xlab = "MAF",
             ylab = "p-value")
})
pmat <- do.call("rbind",lapply(sim_list$sim3,FUN=function(mat){
  mat[,3]
}))
pmat <- matrix(as.vector(pmat),ncol=5)
gg_boxplot(x = rep(seq(0.1,0.5,0.1),each=nrow(pmat)),
           y = as.vector(pmat),
           title = "p-values of permuted traits",
           xlab = "MAF",
           ylab = "p-value")

#sim4
p_list <- lapply(1:length(sim_list$sim1),FUN=function(i){
  plot_mat <- cbind(sim_list$sim1[[i]][,2],
                    sim_list$sim4[[i]][,2])
  max_p <- max(-log10(plot_mat))
  hits <- which(sim_list$sim1[[i]][,1]!=0)
  colvec <- ifelse(1:nrow(plot_mat)%in%hits,"red","black")
  p <- gg_dotplot(x = -log10(plot_mat[,1]),
             y = -log10(plot_mat[,2]),
             ylim=c(0,max_p),
             xlim=c(0,max_p),
             col=colvec,
             xlab = "-log10(p), MB=1",
             ylab = "-log10(p), MB=15",
             main = paste0("simulation ",i))
  add_line(p,slope = 1,offset = 0,color = "blue",linewidth = 0.5)
})
pvec_perm <- unlist(lapply(sim_list$sim4,FUN=function(plot_mat){
  return(plot_mat[plot_mat[,1]==0,3])
}))
gg_histogram(x = pvec_unperm,
             breaks = 20,
             color = "darkgreen",
             xlab = "p-value",
             ylab = "frequency",
             title = "permuted traits, markers without simulated effects")

#sim5
p_list <- lapply(1:length(sim_list$sim5),FUN=function(i){
  plot_mat <- sim_list$sim1[[i]]
  max_p <- max(-log10(plot_mat[,2:3]))
  hits <- which(plot_mat[,1]!=0)
  colvec <- ifelse(1:nrow(plot_mat)%in%hits,"red","black")
  gg_dotplot(x = -log10(plot_mat[,2]),
             y = -log10(plot_mat[,3]),
             ylim=c(0,max_p),
             xlim=c(0,max_p),
             col=colvec,
             xlab = "-log10(p), non permuted",
             ylab = "-log10(p), permuted",
             main = paste0("simulation ",i))
  
})