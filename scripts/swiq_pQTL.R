#####set global parameters######################################################
setwd("/cellfile/cellnet/RF_QTL_package/")
plot_path <- "plots/swiq_exploration/"
trait_path <- "data/mapped_p_traits/"
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
load("data/prot_mapping_data.RData")

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
  
  out <- data.frame(real=p_real)
  saveRDS(object = out,file = paste0(trait_path,"trait_",trait,".rds"))
})
)