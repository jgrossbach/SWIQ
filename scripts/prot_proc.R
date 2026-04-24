library(readr)
load("data/eQtlMappingData160728.RData")
all_trp = read_rds("/cellfile/cellnet/LiPQTL/Data/RMBYcross/all_trp_merged_filtered_new_prot_quant.rds")
prot_samples <- all_trp$trp_prot_merged_maxlfq
rownames(prot_samples) <- prot_samples$Protein
prot_samples <- prot_samples[,-1]

cnames <- gsub(pattern = "__.*",replacement = "",colnames(prot_samples))
strains <- unique(cnames)

phenotype <- sapply(strains,FUN=function(s){
  rowMeans(prot_samples[,cnames==s,drop=F],na.rm=T)
})

number_first <- grepl(pattern = "^[0-9]",x = strains)
strains[number_first] <- paste0("X",strains[number_first])
strains <- gsub("_",".",strains,fixed=T)
colnames(phenotype) <- strains
phenotype <- phenotype[,rownames(genotype)]
phenotype <- t(apply(phenotype,1,FUN=function(x){
  x <- x-mean(x,na.rm=T)
  x <- x/sd(x,na.rm=T)
  return(x)
}))


identical(sort(strains),sort(rownames(genotype)))

save(genotype,phenotype,file = "data/prot_mapping_data.RData")
