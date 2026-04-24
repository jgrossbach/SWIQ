#####set global parameters######################################################
setwd("/cellfile/cellnet/RF_QTL_package/")
plot_path <- "plots/swiq_exploration/"
trait_path <- "data/mapped_albert_traits_small_sample/"
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
folds <- readRDS(file = "data/albert_folds.rds")
folds <- folds$folds
colnames(phenotype) <- rownames(genotype)

#####Map traits#################################################################
done_traits <- list.files(trait_path)
open_traits <- which(!paste0("trait_",1:nrow(phenotype),".rds")%in%
                       done_traits)

invisible(mclapply(X = open_traits,mc.cores = 10,FUN = function(trait){
  print(trait)
  set.seed(trait)
  s <- sample(1:ncol(phenotype),size = 200,replace = F)
  p <- swiq(y = phenotype[trait,s],
            x = genotype[s,],
            num.trees = 20000,
            num.threads = 1,
            min.bucket = 15,
            imp.blocks = 1,
            max.trees = Inf,
            min.trees = 1,
            pick.threshold = 1)
  
  saveRDS(object = data.frame(p=p),file = paste0(trait_path,"trait_",trait,".rds"))
})
)