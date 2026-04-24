#####set global parameters######################################################
setwd("/cellfile/cellnet/RF_QTL_package/")
plot_path <- "plots/lipQTL/"
overwrite <- F

#####load code##################################################################
source("../phosphoQTL/lib/QTL_draft.R") #has to be loaded before rfqtl
library(RFQTL)
library(ranger)
library(RColorBrewer)
library(biomaRt)
source("scripts/swiq.R")
source("scripts/help_fun.R")
source("/cellfile/cellnet/phosphoQTL/lib/redmineTbls.R")

#####load data##################################################################
raw_FT <- readRDS("/cellfile/cellnet/LiPQTL/Data/RMBYcross/LiPAnalyzeR4QTL_20250908/Data/LiPAnalyzeR_LiPPepFT_09022026.rds")
raw_HT <- readRDS("/cellfile/cellnet/LiPQTL/Data/RMBYcross/LiPAnalyzeR4QTL_20250908/Data/LiPAnalyzeR_LiPPepHT_09022026.rds")
load("data/lip_mapping_data.RData")
all_lip_merged <- readRDS("/cellfile/cellnet/LiPQTL/Data/RMBYcross/all_lip_merged_filtered.rds")
pep2prot <- rbind(all_lip_merged$lip_pep_merged_ft[,1:2],all_lip_merged$lip_pep_merged_ht[,1:2])
pep2prot <- as.data.frame(pep2prot)
pep2prot <- pep2prot[pep2prot[,1]%in%rownames(phenotype),]
rownames(pep2prot) <- pep2prot$Peptide
pep2prot <- pep2prot[rownames(phenotype),]
pep2prot <- pep2prot$Protein
names(pep2prot) <- rownames(phenotype)

load("data/plotPack.RData")
names(chr_length) <- gsub(pattern = "chr",
                          replacement = "",
                          x = names(chr_length))
chrLen <- chr_length
chrCol <- "grey"

chrLenCum <- sapply(1:length(chrLen),FUN=function(x){
  sum(chrLen[1:x])
})

geno <- read.table("data/genotype_for_mapping.tsv",header=T)
chrVec <- geno[,1]
mloc <- geno[,1:3]
geno <- geno[,4:ncol(geno)]
geno <- t(geno)

#biomart
mart <- useMart("ensembl")
ensembl <- useDataset("scerevisiae_gene_ensembl", mart = mart)

#peptide locations
pep_loc <- readRDS("/cellfile/cellnet/LiPQTL/Data/RMBYcross/LiPAnalyzeR4QTL_20250908/Data/LiPAnalyzeR_annotPepProt_cleaned_18022026.rds")
rownames(pep_loc) <- pep_loc$Peptide
#####results####################################################################
trait_path <- "data/mapped_lip_traits/"
files <- list.files(trait_path)
p <- t(sapply(1:length(files),FUN=function(i){
  fname <- paste0(trait_path,"trait_",i,".rds")
  out <- readRDS(fname)
  return(out[,1])
}))
rownames(p) <- rownames(phenotype)
q_HT <- p[rownames(phenotype)%in%rownames(raw_HT$modelResid),]
q_HT[] <- p.adjust(q_HT,method="fdr")

q_FT <- p[rownames(phenotype)%in%rownames(raw_FT$modelResid),]
q_FT[] <- p.adjust(q_FT,method="fdr")

fdrs <- c(1,5,10,15,20,25) 
qtl_HT <- lapply(fdrs,FUN=function(fdr){
  QTLgrouper(pmat = q_HT,
             sigThreshold = fdr/100,
             corThreshold = 0.8,
             distThreshold = 9,
             genotype = geno,
             chrVec = chrVec)
})
names(qtl_HT) <- paste0("FDR",fdrs)

qtl_FT <- lapply(fdrs,FUN=function(fdr){
  QTLgrouper(pmat = q_FT,
             sigThreshold = fdr/100,
             corThreshold = 0.8,
             distThreshold = 9,
             genotype = geno,
             chrVec = chrVec)
})
names(qtl_FT) <- paste0("FDR",fdrs)


#call local qtl
# get gene prosition (chr coordinates)
gff <- read.table(file = "/cellfile/cellnet/phosphoQTL/Saccharomyces_cerevisiae/RM/genes_RM_CDSonly_refined.gtf",stringsAsFactors = F)
g <- read.table("data/genotype_for_mapping.tsv",header=T,check.names = F)
colnames(g)[1]<-"chr"
CDS_position<-t(sapply(unique(gff[,9]),function(gene){
  sel<-gff[,9]==gene
  pos<-c(gff[sel,4],gff[sel,5])
  pos<-range(pos)
  chr<-gff[sel,1][1]
  return(list(chr,pos[1],pos[2]))
}))
gene<-CDS_position
gene<-data.frame(chr=unlist(gene[,1]),
                 start=unlist(gene[,2]),
                 end=unlist(gene[,3]),
                 stringsAsFactors = F)
gene$med <- rowMeans(cbind(gene$start,gene$end))

# get gene position (relative to genetic markers)
geneMarker<-sapply(1:nrow(gene),function(i){
  before <- max(which(as.character(g$chr)==as.character(gene$chr[i]) &
                        g$end<gene$med[i]))
  if(is.infinite(before)) {before =min(which(as.character(g$chr)==as.character(gene$chr[i])))}
  after <- min(which(as.character(g$chr)==as.character(gene$chr[i]) &
                       g$start>gene$med[i]))
  if(is.infinite(after)) {after =max(which(as.character(g$chr)==as.character(gene$chr[i])))}
  return(c(before,after))
})
gene <- cbind(gene,marker_1=geneMarker[1,],marker_2=geneMarker[2,],stringsAsFactors=F)
cor.all <- cor(genotype,use="pair")
qtl_FT <-lapply(qtl_FT,FUN=function(qtl_list){
  lapply(qtl_list,function(qtl){
    target <- rownames(q_FT)[qtl$target]
    qtl$pep <- target
    target <- pep2prot[target]
    qtl$prot <- target
    M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
    if(any(is.infinite(M))){
      qtl$local <- F
    }else{
      COR<-max(cor.all[M[1]:M[2],
                       unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
      qtl$local <- COR>0.8
    }
    return(qtl)
  })
})
qtl_HT <-lapply(qtl_HT,FUN=function(qtl_list){
  lapply(qtl_list,function(qtl){
    target <- rownames(q_HT)[qtl$target]
    qtl$pep <- target
    target <- pep2prot[target]
    qtl$prot <- target
    M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
    if(any(is.infinite(M))){
      qtl$local <- F
    }else{
      COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
      qtl$local <- COR>0.8
    }
    return(qtl)
  })
})

#summarize the number of qtls
qtl_tables <- lapply(list(qtl_FT,qtl_HT),FUN=function(qtl){
  
  qtl_numbers <- data.frame(FDR=fdrs)
  qtl_numbers$n_QTL <- sapply(qtl,length)
  qtl_numbers$n_peps <- sapply(qtl,FUN=function(qtl_list){
    length(unique(sapply(qtl_list,FUN=function(qtl){
      qtl$pep
    })))
  })
  qtl_numbers$n_proteins <- sapply(qtl,FUN=function(qtl_list){
    length(unique(sapply(qtl_list,FUN=function(qtl){
      qtl$prot
    })))
  })
  qtl_numbers$n_QTL_local <- sapply(qtl,FUN=function(qtl_list){
    ct <- 0
    local <- sapply(qtl_list,FUN=function(qtl){
      ct <<- ct+1
      qtl$local
    })
    sum(local)
  })
  qtl_numbers$n_peps_local <- sapply(qtl,FUN=function(qtl_list){
    ct <- 0
    local <- sapply(qtl_list,FUN=function(qtl){
      ct <<- ct+1
      qtl$local
    })
    target <- sapply(qtl_list,FUN=function(qtl){
      qtl$target
    })
    length(unique(target[local]))
  })
  qtl_numbers$n_proteins_local <- sapply(qtl,FUN=function(qtl_list){
    ct <- 0
    local <- sapply(qtl_list,FUN=function(qtl){
      ct <<- ct+1
      qtl$local
    })
    target <- sapply(qtl_list,FUN=function(qtl){
      qtl$prot
    })
    length(unique(target[local]))
  })
  
  rownames(qtl_numbers) <- NULL
  return(qtl_numbers)
})


sQTL_results <- list(q_HT=q_HT,
                     q_FT=q_FT,
                     qtl_HT=qtl_HT,
                     qtl_FT=qtl_FT,
                     qtl_tables=qtl_tables,
                     pep2prot=pep2prot,
                     genotype=genotype,
                     phenotype=phenotype)
saveRDS(file = "data/sQTL_results.rds",object = sQTL_results)

qtl_targets <- lapply(list(qtl_FT$FDR10,qtl_HT$FDR10),FUN=function(qtl_list){
  targets <- data.frame(protein=sapply(qtl_list,FUN=function(qtl){
    qtl$prot
  }),
  local=sapply(qtl_list,FUN=function(qtl){
    qtl$local
  }))
  return(targets)
})
names(qtl_targets) <- c("FT","HT")

uniprots <- unique(c(qtl_targets$FT$protein,qtl_targets$HT$protein))
overlap <- table(uniprots%in%qtl_targets$FT[qtl_targets$FT$local,1],
                 uniprots%in%qtl_targets$HT[qtl_targets$HT$local,1])
fisher.test(overlap)
printRedmine(overlap)

#make tables with all relevant qtl infos
HT_table <- qtl_targets$HT
prot_infos <- getBM(
  attributes = c("ensembl_gene_id","uniprotswissprot","external_gene_name"),
  filters = "ensembl_gene_id",
  values = HT_table$protein,
  mart = ensembl
)
rownames(prot_infos) <- prot_infos$ensembl_gene_id

HT_table <- cbind(HT_table,prot_infos[HT_table$protein,])

qtl_info <- t(sapply(qtl_HT$FDR10,FUN=function(qtl){
  chr <- mloc[qtl$predictors[1,1],1]
  start <- min(mloc[as.vector(qtl$predictors),2:3])
  end <- max(mloc[as.vector(qtl$predictors),2:3])
  return(c(chr,start,end,end-start+1,qtl$pep))
}))
qtl_info <- as.data.frame(qtl_info)
colnames(qtl_info) <- c("QTL_chr","QTL_start","QTL_end","QTL_size","peptide")
HT_table <- cbind(HT_table,qtl_info)
HT_table$peptide_type <- "HT"
pep_infos <- pep_loc[HT_table$peptide,c("startPosition","endPosition")]
colnames(pep_infos) <- c("peptide_start","peptide_end")
HT_table <- cbind(HT_table,pep_infos)
HT_table$protein <- NULL

FT_table <- qtl_targets$FT
prot_infos <- getBM(
  attributes = c("ensembl_gene_id","uniprotswissprot","external_gene_name"),
  filters = "ensembl_gene_id",
  values = FT_table$protein,
  mart = ensembl
)
rownames(prot_infos) <- prot_infos$ensembl_gene_id

FT_table <- cbind(FT_table,prot_infos[FT_table$protein,])
qtl_info <- t(sapply(qtl_FT$FDR10,FUN=function(qtl){
  chr <- mloc[qtl$predictors[1,1],1]
  start <- min(mloc[as.vector(qtl$predictors),2:3])
  end <- max(mloc[as.vector(qtl$predictors),2:3])
  return(c(chr,start,end,end-start+1,qtl$pep))
}))
qtl_info <- as.data.frame(qtl_info)
colnames(qtl_info) <- c("QTL_chr","QTL_start","QTL_end","QTL_size","peptide")
FT_table <- cbind(FT_table,qtl_info)
FT_table$peptide_type <- "FT"
pep_infos <- pep_loc[FT_table$peptide,c("startPosition","endPosition")]
colnames(pep_infos) <- c("peptide_start","peptide_end")
FT_table <- cbind(FT_table,pep_infos)
FT_table$protein <- NULL

sQTL_table <- rbind(HT_table,FT_table)
sQTL_table <- sQTL_table[,c("ensembl_gene_id",
                            "uniprotswissprot",
                            "external_gene_name",
                            "peptide",
                            "peptide_type",
                            "peptide_start",
                            "peptide_end",
                            "QTL_chr",
                            "QTL_start",
                            "QTL_end",
                            "QTL_size",
                            "local")]
write.table(x = sQTL_table,
            file = "data/sQTL.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)

#####plots######################################################################

png(filename = paste0(plot_path,"qtl10_FT.png"),height = 1300,width=1300,res = 200)
QTLplotter(QTLlist = qtl_FT$FDR10,
           targetLocs = geneLocs[pep2prot[rownames(q_FT)]],
           predictorLocs = markerLocs,
           chrLen = chr_length,
           xlab = "QTL position",
           ylab = "target-gene position",
           qtl.lwd=5,
           labcex = 1,
           col=adjustcolor(col = "black",alpha.f = 0.2),
           chrCol = rgb(0,0,0,0.3),
           chrLty=1,
           main="sQTL for FT peptides")
dev.off()

png(filename = paste0(plot_path,"qtl10_HT.png"),height = 1300,width=1300,res = 200)
QTLplotter(QTLlist = qtl_HT$FDR10,
           targetLocs = geneLocs[pep2prot[rownames(q_HT)]],
           predictorLocs = markerLocs,
           chrLen = chr_length,
           xlab = "QTL position",
           ylab = "target-gene position",
           qtl.lwd=5,
           labcex = 1,
           col=adjustcolor(col = "black",alpha.f = 0.2),
           chrCol = rgb(0,0,0,0.3),
           chrLty=1,
           main="sQTL for HT peptides")
dev.off()

# qtl_lead <- lapply(qtl$FDR10,FUN=function(qtl){
#   qtl$predictors <- matrix(rep(qtl$mostSignificantPredictor,2),nrow=1,ncol=2)
#   colnames(qtl$predictors) <- c("start","end")
#   return(qtl)
# })
# 
# QTLplotter(QTLlist = qtl_lead,
#            targetLocs = geneLocs[pep2prot[rownames(q)]],
#            predictorLocs = markerLocs,
#            chrLen = chr_length,
#            xlab = "QTL position",
#            ylab = "target-gene position",
#            qtl.lwd=5,
#            labcex = 1,
#            col=adjustcolor(col = "black",alpha.f = 0.2),
#            chrCol = rgb(0,0,0,0.3),
#            chrLty=1,
#            main="sQTL")

png(filename = paste0(plot_path,"qtl10_FT_profile.png"),height = 800,width=1300,res = 200)
plot(y = colSums(q_FT<0.1),
     x = markerLocs,
     type="h",
     xaxt="n",
     ylab="number of peptides",
     xlab="QTL location",
     main="FT peptides")
for(i in 1:(length(chrLenCum)-1)){
  abline(v=chrLenCum[i],lty=2,col="grey")
}
borders <- c(0,chrLenCum)
for(i in 1:length(chrLen)){
  mtext(names(chrLen[i]),
        side=1,
        at=mean(borders[i:(i+1)]),
        cex=par()$cex*1*0.75)
}
dev.off()

png(filename = paste0(plot_path,"qtl10_HT_profile.png"),height = 800,width=1300,res = 200)
plot(y = colSums(q_HT<0.1),
     x = markerLocs,
     type="h",
     xaxt="n",
     ylab="number of peptides",
     xlab="QTL location",
     main="HT peptides")
for(i in 1:(length(chrLenCum)-1)){
  abline(v=chrLenCum[i],lty=2,col="grey")
}
borders <- c(0,chrLenCum)
for(i in 1:length(chrLen)){
  mtext(names(chrLen[i]),
        side=1,
        at=mean(borders[i:(i+1)]),
        cex=par()$cex*1*0.75)
}
dev.off()


sapply(qtl_tables,printRedmine)
