#####set global parameters######################################################
setwd("/cellfile/cellnet/RF_QTL_package/")
plot_path <- "plots/swiq_exploration/"
trait_path <- "data/mapped_e_traits/"
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

#####Map traits#################################################################
done_traits <- list.files(trait_path)
open_traits <- which(!paste0("trait_",1:nrow(phenotype),".rds")%in%
                       done_traits)
for(trait in open_traits){
  p_real <- swiq(y = phenotype[trait,],
                 x = genotype[,1:3593],
                 num.trees = 2500,
                 num.threads = 15,
                 min.bucket = 15,
                 imp.blocks = 4)
  set.seed(trait)
  idx <- sample(1:nrow(genotype),size=nrow(genotype))
  p_perm <- swiq(y = phenotype[trait,idx],
                 x = genotype[,1:3593],
                 num.trees = 2500,
                 num.threads = 15,
                 min.bucket = 15,
                 imp.blocks = 4)
  out <- data.frame(real=p_real,
                    perm=p_perm)
  saveRDS(object = out,file = paste0(trait_path,"trait_",trait,".rds"))
}