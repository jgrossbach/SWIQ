load("data/eQtlMappingData160728.RData")
#raw <- readRDS("data/LiPAnalyzeR_LiPPepFT_27012026.rds")

# FT: "/cellfile/cellnet/LiPQTL/Data/RMBYcross/LiPAnalyzeR4QTL_20250908/Data/LiPAnalyzeR_LiPPepFT_09022026.rds"
# HT: "/cellfile/cellnet/LiPQTL/Data/RMBYcross/LiPAnalyzeR4QTL_20250908/Data/LiPAnalyzeR_LiPPepHT_09022026.rds"

raw_FT <- readRDS("/cellfile/cellnet/LiPQTL/Data/RMBYcross/LiPAnalyzeR4QTL_20250908/Data/LiPAnalyzeR_LiPPepFT_09022026.rds")
raw_HT <- readRDS("/cellfile/cellnet/LiPQTL/Data/RMBYcross/LiPAnalyzeR4QTL_20250908/Data/LiPAnalyzeR_LiPPepHT_09022026.rds")

identical(colnames(raw_FT$modelResid),colnames(raw_HT$modelResid))

cnames <- gsub(pattern = "__.*",replacement = "",colnames(raw_FT$modelResid))
strains <- unique(cnames)
phenotype_FT <- sapply(strains,FUN=function(s){
  rowMeans(raw_FT$modelResid[,cnames==s,drop=F],na.rm=T)
})

phenotype_HT <- sapply(strains,FUN=function(s){
  rowMeans(raw_HT$modelResid[,cnames==s,drop=F],na.rm=T)
})

phenotype <- rbind(phenotype_FT,phenotype_HT)

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

save(genotype,phenotype,file = "data/lip_mapping_data.RData")
