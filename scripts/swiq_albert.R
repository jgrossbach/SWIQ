#####set global parameters######################################################
setwd("/cellfile/cellnet/RF_QTL_package/")
plot_path <- "plots/swiq_exploration/"
trait_path <- "data/mapped_albert_traits/"
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
mapping_data <- readRDS(file = "data/albert_mapping_data.rds")
genotype <- mapping_data$genotype
phenotype <- mapping_data$phenotype

#####Map traits#################################################################
done_traits <- list.files(trait_path)
open_traits <- which(!paste0("trait_",1:nrow(phenotype),".rds")%in%
                       done_traits)

invisible(mclapply(X = open_traits,mc.cores = 10,FUN = function(trait){
  p_real <- swiq(y = phenotype[trait,],
                 x = genotype,
                 num.trees = 20000,
                 num.threads = 1,
                 min.bucket = 15,
                 imp.blocks = 1,
                 max.trees = Inf,
                 min.trees = 1,
                 pick.threshold = 1)
  y_perm <- phenotype[trait,]
  y_perm <- sample(y_perm,size=length(y_perm))
  p_perm <- swiq(y = y_perm,
                 x = genotype,
                 num.trees = 20000,
                 num.threads = 1,
                 min.bucket = 15,
                 imp.blocks = 1,
                 max.trees = Inf,
                 min.trees = 1,
                 pick.threshold = 1)
  
  out <- data.frame(real=p_real,perm=p_perm)
  saveRDS(object = out,file = paste0(trait_path,"trait_",trait,".rds"))
})
)