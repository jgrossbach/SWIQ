mapping_data <- readRDS("data/albert_mapping_data.rds")
genotype <- mapping_data$genotype
phenotype <- mapping_data$phenotype
colnames(phenotype) <- rownames(genotype)
set.seed(3232)

strain_fold <- rep(1:2,length.out=nrow(genotype))
strain_fold <- strain_fold[sample(1:nrow(genotype))]
fcs <- lapply(1:max(strain_fold),FUN=function(i){
  i <<- i
  strains <- rownames(genotype)[which(strain_fold==i)]
  phenotype <- phenotype-rowMeans(phenotype,na.rm=T)
  phenotype <- phenotype/apply(phenotype,1,sd,na.rm=T)
  fc_by_m <- sapply(1:ncol(genotype),FUN=function(m){
    m <<- m
    s0 <- strains[which(genotype[strains,m]==0)]
    s1 <- strains[which(genotype[strains,m]==1)]
    fc <- rowMeans(phenotype[,s1],na.rm=T)-rowMeans(phenotype[,s0],na.rm=T)
    return(fc)
  })
  return(fc_by_m)
})
fold_info <- list(folds=strain_fold,fcs=fcs)
saveRDS(object = fold_info,file = "data/albert_folds.rds")
