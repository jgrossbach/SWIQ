#####set global parameters######################################################
setwd("/cellfile/cellnet/RF_QTL_package/")
plot_path <- "plots/swiq_exploration/"
trait_path <- "data/mapped_e_traits_skimmed/"
overwrite <- F

#####load code##################################################################
source("../phosphoQTL/lib/QTL_draft.R") #has to be loaded before rfqtl
library(RFQTL)
library(ranger)
library(RColorBrewer)
library(parallel)
library(abind)
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

#####Map traits#################################################################
done_traits <- list.files(trait_path)
open_traits <- which(!paste0("trait_",1:nrow(phenotype),".rds")%in%
                       done_traits)

invisible(mclapply(X = open_traits,mc.cores = 15,FUN = function(trait){
  p_real <- swiq(y = phenotype[trait,],
                 x = genotype[,1:3593],
                 num.trees = 2500,
                 num.threads = 1,
                 min.bucket = 15,
                 imp.blocks = 4,
                 pick.threshold = 0.001)
  
  p_str <- swiq(y = phenotype[trait,],
                x = genotype[,1:3593],
                num.trees = 2500,
                num.threads = 1,
                min.bucket = 15,
                imp.blocks = 4,
                pick.threshold = 0.001,
                pop.str = as.data.frame(genotype[,3594:3600]))
  
  out <- data.frame(real=p_real,
                    str=p_str)
  saveRDS(object = out,file = paste0(trait_path,"trait_",trait,".rds"))
})
)