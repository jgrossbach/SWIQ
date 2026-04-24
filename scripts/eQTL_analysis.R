#####set global parameters######################################################
setwd("/cellfile/cellnet/RF_QTL_package/")
plot_path <- "plots/swiq_exploration/"
overwrite <- F

#####load code##################################################################
source("../phosphoQTL/lib/QTL_draft.R") #has to be loaded before rfqtl
library(RFQTL)
library(ranger)
library(RColorBrewer)
library(biomaRt)
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
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]
names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")

geno <- read.table("data/genotype_for_mapping.tsv",header=T)
chrVec <- geno[,1]
mloc <- geno[,1:3]
geno <- geno[,4:ncol(geno)]
geno <- t(geno)

#biomart
mart <- useMart("ensembl")
ensembl <- useDataset("scerevisiae_gene_ensembl", mart = mart)

#####results####################################################################
eQTL_path <- "data/mapped_e_traits_skimmed/"
files <- list.files(eQTL_path)
eQTL_p <- t(sapply(1:length(files),FUN=function(i){
  fname <- paste0(eQTL_path,"trait_",i,".rds")
  out <- readRDS(fname)
  return(out[,1])
}))
eQTL_q <- eQTL_p
eQTL_q[] <- p.adjust(eQTL_p,method="fdr")
rownames(eQTL_q) <- rownames(phenotype)

fdrs <- c(1,5,10,15,20,25) 
eQTL <- lapply(fdrs,FUN=function(fdr){
  QTLgrouper(pmat = eQTL_q,
             sigThreshold = fdr/100,
             corThreshold = 0.8,
             distThreshold = 9,
             genotype = geno,
             chrVec = chrVec)
})
names(eQTL) <- paste0("FDR",fdrs)

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
eQTL <-lapply(eQTL,FUN=function(qtl_list){
  lapply(qtl_list,function(qtl){
    target <- rownames(eQTL_q)[qtl$target]
    #qtl$pep <- target
    #target <- pep2prot[target]
    qtl$gene <- target
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

targets <- data.frame(protein=sapply(eQTL$FDR10,FUN=function(qtl){
  qtl$gene
}),
local=sapply(eQTL$FDR10,FUN=function(qtl){
  qtl$local
}))
qtl_table <- targets
prot_infos <- getBM(
  attributes = c("ensembl_gene_id","uniprotswissprot","external_gene_name"),
  filters = "ensembl_gene_id",
  values = rownames(eQTL_q),
  mart = ensembl
)
prot_infos <- lapply(unique(qtl_table$protein),FUN=function(gene){
  prot_infos[which(prot_infos$ensembl_gene_id==gene)[1],]
})
prot_infos <- do.call("rbind",prot_infos)
prot_infos <- prot_infos[!is.na(prot_infos$ensembl_gene_id),]
rownames(prot_infos) <- prot_infos$ensembl_gene_id

qtl_table <- cbind(qtl_table,prot_infos[qtl_table$protein,])

qtl_info <- t(sapply(eQTL$FDR10,FUN=function(qtl){
  chr <- mloc[qtl$predictors[1,1],1]
  start <- min(mloc[as.vector(qtl$predictors),2:3])
  end <- max(mloc[as.vector(qtl$predictors),2:3])
  return(c(chr,start,end,end-start+1))
}))
qtl_info <- as.data.frame(qtl_info)
colnames(qtl_info) <- c("QTL_chr","QTL_start","QTL_end","QTL_size")
qtl_table <- cbind(qtl_table,qtl_info)
qtl_table$protein <- NULL

eQTL_table <- qtl_table[,c("ensembl_gene_id",
                            "uniprotswissprot",
                            "external_gene_name",
                            "QTL_chr",
                            "QTL_start",
                            "QTL_end",
                            "QTL_size",
                            "local")]
eQTL_table <- eQTL_table[!is.na(eQTL_table$ensembl_gene_id),]
write.table(x = eQTL_table,
            file = "data/eQTL.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)

local_effects <- sapply(rownames(phenotype),FUN=function(g){
  m <- ceiling(mean(as.numeric(gene[g,5:6])))
  p <- scale(phenotype[g,],center = T,scale = T)[,1]
  effect <- mean(p[which(genotype[,m]==1)],na.rm=T)-mean(p[which(genotype[,m]==0)],na.rm=T)
  return(effect)
})

eQTL_results <- list(q=eQTL_q,
                     qtl=eQTL,
                     genotype=genotype,
                     phenotype=phenotype,
                     local_effects=local_effects)
saveRDS(file = "data/eQTL_results.rds",object = eQTL_results)



#####plots######################################################################
png(filename = paste0(plot_path,"eQTL_map.png"),height = 1300,width=1300,res = 200)
QTLplotter(QTLlist = eQTL$FDR10,
           targetLocs = geneLocs[rownames(qv)],
           predictorLocs = markerLocs,
           chrLen = chr_length,
           xlab = "eQTL position",
           ylab = "target-gene position",
           qtl.lwd=5,
           labcex = 1,
           col=adjustcolor(col = colors["eQTL"],alpha.f = 0.2),
           chrCol = rgb(0,0,0,0.3),
           chrLty=1,
           main="eQTL from phosphoQTL data")
dev.off()


qtl_lead <- lapply(eQTL$FDR10,FUN=function(qtl){
  qtl$predictors <- data.frame(start=qtl$mostSignificantPredictor,
                               end=qtl$mostSignificantPredictor)
  return(qtl)
})
png(filename = paste0(plot_path,"eQTL_phosphoQTL_map_lead.png"),height = 1300,width=1300,res = 200)
QTLplotter(QTLlist = qtl_lead,
           targetLocs = geneLocs[rownames(qv)],
           predictorLocs = markerLocs,
           chrLen = chr_length,
           xlab = "eQTL position",
           ylab = "target-gene position",
           qtl.lwd=5,
           labcex = 1,
           col=adjustcolor(col = colors["eQTL"],alpha.f = 0.2),
           chrCol = rgb(0,0,0,0.3),
           chrLty=1,
           main="SWIQ")
dev.off()

qtl_lead <- lapply(QtlList$FDR10,FUN=function(qtl){
  qtl$predictors <- data.frame(start=qtl$mostSignificantPredictor,
                               end=qtl$mostSignificantPredictor)
  return(qtl)
})
png(filename = paste0(plot_path,"eQTL_phosphoQTL_map_lead_legacy.png"),height = 1300,width=1300,res = 200)
QTLplotter(QTLlist = qtl_lead,
           targetLocs = geneLocs[rownames(qv)],
           predictorLocs = markerLocs,
           chrLen = chr_length,
           xlab = "eQTL position",
           ylab = "target-gene position",
           qtl.lwd=5,
           labcex = 1,
           col=adjustcolor(col = colors["eQTL"],alpha.f = 0.2),
           chrCol = rgb(0,0,0,0.3),
           chrLty=1,
           main="legacy method")
dev.off()

chrLenCum <- sapply(1:length(chr_length),FUN=function(x){
  sum(chr_length[1:x])
})
png(filename = paste0(plot_path,"eQTL_profile.png"),height = 800,width=1300,res = 200)
plot(y = colSums(eQTL_q<0.1),
     x = markerLocs,
     type="h",
     xaxt="n",
     ylab="number of traits",
     xlab="QTL location",
     main="transcripts")
for(i in 1:(length(chrLenCum)-1)){
  abline(v=chrLenCum[i],lty=2,col="grey")
}
borders <- c(0,chrLenCum)
for(i in 1:length(chr_length)){
  mtext(names(chr_length[i]),
        side=1,
        at=mean(borders[i:(i+1)]),
        cex=par()$cex*1*0.75)
}
dev.off()


eQTL_lead <- lapply(eQTL$FDR10,FUN=function(qtl){
  qtl$predictors <- matrix(rep(qtl$mostSignificantPredictor,2),nrow=1,ncol=2)
  colnames(qtl$predictors) <- c("start","end")
  return(qtl)
})

scaled_q <- sapply(eQTL_lead,FUN=function(qtl){
  eQTL_q[qtl$target,qtl$mostSignificantPredictor]
})
scaled_q <- -log10(scaled_q)
scaled_q <- scaled_q/(max(scaled_q)/2)
scaled_q[scaled_q>1] <- 1
scaled_q <- sapply(scaled_q,FUN=function(x){
  adjustcolor(col = colors["eQTL"],alpha.f = x)
})
QTLplotter(QTLlist = eQTL_lead,
           targetLocs = geneLocs[rownames(qv)],
           predictorLocs = markerLocs,
           chrLen = chr_length,
           xlab = "eQTL position",
           ylab = "target-gene position",
           qtl.lwd=5,
           labcex = 1,
           col=adjustcolor(col = colors["eQTL"],alpha.f = 0.2),
           chrCol = rgb(0,0,0,0.3),
           chrLty=1,
           main="swiq results")

table(table(sapply(eQTL$FDR10,FUN=function(qtl){qtl$target})))
table(table(sapply(QtlList$FDR10,FUN=function(qtl){qtl$target})))


QTLplotter(QTLlist = QtlList$FDR10,
           targetLocs = geneLocs[rownames(qv)],
           predictorLocs = markerLocs,
           chrLen = chr_length,
           xlab = "eQTL position",
           ylab = "target-gene position",
           qtl.lwd=5,
           labcex = 1,
           col=adjustcolor(col = colors["eQTL"],alpha.f = 0.2),
           chrCol = rgb(0,0,0,0.3),
           chrLty=1,
           main="old results")
