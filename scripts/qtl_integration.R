#####set global parameters######################################################
plot_path <- "plots/phosphoQTL_benchmark/"

#####load code##################################################################
source("/cellfile/cellnet/phosphoQTL/lib/redmineTbls.R")
source("scripts/help_fun.R")
library(readxl)
library(ggplot2)
library(gridExtra)

#####load data##################################################################
eQTL <- read.table("data/eQTL.tsv",sep="\t",header=T)
pQTL <- read.table("data/pQTL.tsv",sep="\t",header=T)
pQTL2 <- read.table("data/pQTL_phosphoQTL.tsv",sep="\t",header=T)
ptQTL <- read.table("data/ptQTL.tsv",sep="\t",header=T)
sQTL <- read.table("data/sQTL.tsv",sep="\t",header=T)
phrQTL <- read.table("data/phrQTL.tsv",sep="\t",header=T)
legacy_results <- readRDS("data/phosphoQTL_results/qtl_tables.rds")
legacy_results <- legacy_results[1:4]
swiq_results <- list(eQTL=eQTL,ptQTL=ptQTL,pQTL=pQTL2,phrQTL=phrQTL)
alberts_eQTL <- as.data.frame(read_excel(path = "data/albert2018/SI_Data_04_eQTL.xlsx"))
alberts_local <- as.data.frame(read_excel(path = "data/albert2018/SI_Data_07_ASE_localeQTL.xlsx"))


#####analysis###################################################################
layers <- c("eQTL","pQTL","pQTL2","ptQTL")
n_traits <- c(5429,3856,1862,1857)
names(n_traits) <- layers
qtl_list <- list(eQTL=eQTL,
                 pQTL=pQTL,
                 pQTL2=pQTL2,
                 ptQTL=ptQTL)

qtl_table <- sapply(layers,FUN=function(qtl_type){
  n <- n_traits[qtl_type]
  qtls <- qtl_list[[qtl_type]]
  n_qtl <- nrow(qtls)
  n_local <- sum(qtls$local)
  n_distant <- n_qtl-n_local
  n_targets <- length(unique(qtls$ensembl_gene_id))
  qtl_per_target <- n_qtl/n_targets
  qtl_per_trait <- n_qtl/n
  return(c(n_qtl,
           n_local,
           n_distant,
           n_targets,
           n_local/n,
           n,
           qtl_per_target,
           qtl_per_trait))
})
rownames(qtl_table) <- NULL
qtl_table <- data.frame(metric=c("number of QTL",
                                 "number of local QTL",
                                 "number of distant QTL",
                                 "number of traits with a QTL",
                                 "proportion of traits with a local QTL",
                                 "total number of traits",
                                 "number of QTL per trait with at least one QTL",
                                 "QTL per trait"),
                        qtl_table)

printRedmine(qtl_table)

layer_comparisons <- expand.grid(layers,layers)
layer_comparisons <- layer_comparisons[layer_comparisons[,1]!=layer_comparisons[,2],]
layer_comparisons <- t(apply(layer_comparisons,1,sort))
layer_comparisons <- unique(layer_comparisons)
rownames(layer_comparisons) <- NULL
cont_tables <- lapply(1:nrow(layer_comparisons),FUN=function(i){
  l1 <- qtl_list[[layer_comparisons[i,1]]]
  l2 <- qtl_list[[layer_comparisons[i,2]]]
  # genes <- unique(c(l1$ensembl_gene_id,
  #                   l2$ensembl_gene_id))
  genes <- intersect(l1$ensembl_gene_id,l2$ensembl_gene_id)
  l1_local <- genes%in%l1[l1$local,"ensembl_gene_id"]
  l2_local <- genes%in%l2[l2$local,"ensembl_gene_id"]
  out <- table(l1_local,l2_local)
  return(out)
})
names(cont_tables) <- apply(layer_comparisons,1,paste0,collapse="")

enrichment <- t(sapply(cont_tables,FUN=function(mat){
  p <- fisher.test(mat,alternative="greater")$p.value
  s1 <- mat[2,2]/sum(mat[,2])
  s2 <- mat[2,2]/sum(mat[2,])
  return(c(p,s1,s2))
}))

enrichment <- data.frame(layer_comparisons[,2:1],enrichment)
colnames(enrichment) <- c("A","B","p","A in B","B in A")
rownames(enrichment) <- NULL
enrichment$A <- as.character(enrichment$A)
enrichment$B <- as.character(enrichment$B)
enrichment[enrichment=="pQTL"] <- "pQTL (sQTL)"
enrichment[enrichment=="pQTL2"] <- "pQTL (phosphoQTL)"
enrichment <- enrichment[order(enrichment$A,enrichment$B),]
enrichment[,4:5] <- apply(enrichment[,4:5],2,round,digits=2)
enrichment[,3] <- formatC(enrichment[,3], format = "e", digits = 1)

printRedmine(enrichment)

window_padding <- 10000
qtl_overlap <- t(sapply(1:nrow(layer_comparisons),FUN=function(i){
  l1 <- qtl_list[[layer_comparisons[i,1]]]
  l2 <- qtl_list[[layer_comparisons[i,2]]]
  genes <- intersect(l1$ensembl_gene_id,l2$ensembl_gene_id)
  shared_a_b <- sapply(genes,FUN=function(g){
    l1s <- l1[l1$ensembl_gene_id==g,,drop=F]
    l2s <- l2[l2$ensembl_gene_id==g,,drop=F]
    same_chr <- outer(l1s$QTL_chr,l2s$QTL_chr,"==")
    not_early <- outer(l1s$QTL_end+window_padding,l2s$QTL_start-window_padding,">=")
    not_late <- outer(l1s$QTL_start-window_padding,l2s$QTL_end+window_padding,"<=")
    overlap <- same_chr&not_late&not_early
    overlap <- mean(rowSums(overlap)>0)
    return(overlap)
  })
  shared_b_a <- sapply(genes,FUN=function(g){
    l1s <- l1[l1$ensembl_gene_id==g,,drop=F]
    l2s <- l2[l2$ensembl_gene_id==g,,drop=F]
    same_chr <- outer(l2s$QTL_chr,l1s$QTL_chr,"==")
    not_early <- outer(l2s$QTL_end+window_padding,l1s$QTL_start-window_padding,">=")
    not_late <- outer(l2s$QTL_start-window_padding,l1s$QTL_end+window_padding,"<=")
    overlap <- same_chr&not_late&not_early
    overlap <- mean(rowSums(overlap)>0)
    return(overlap)
  })
  return(c(mean(shared_a_b),mean(shared_b_a)))
}))
qtl_overlap_perm <- t(sapply(1:nrow(layer_comparisons),FUN=function(i){
  l1 <- qtl_list[[layer_comparisons[i,1]]]
  l2 <- qtl_list[[layer_comparisons[i,2]]]
  genes <- intersect(l1$ensembl_gene_id,l2$ensembl_gene_id)
  shared_a_b <- sapply(genes,FUN=function(g){
    l1s <- l1[l1$ensembl_gene_id==g,,drop=F]
    g2 <- sample(genes,size = 1)
    l2s <- l2[l2$ensembl_gene_id==g2,,drop=F]
    same_chr <- outer(l1s$QTL_chr,l2s$QTL_chr,"==")
    not_early <- outer(l1s$QTL_end+window_padding,l2s$QTL_start-window_padding,">=")
    not_late <- outer(l1s$QTL_start-window_padding,l2s$QTL_end+window_padding,"<=")
    overlap <- same_chr&not_late&not_early
    overlap <- mean(rowSums(overlap)>0)
    return(overlap)
  })
  shared_b_a <- sapply(genes,FUN=function(g){
    l1s <- l1[l1$ensembl_gene_id==g,,drop=F]
    g2 <- sample(genes,size = 1)
    l2s <- l2[l2$ensembl_gene_id==g2,,drop=F]
    same_chr <- outer(l2s$QTL_chr,l1s$QTL_chr,"==")
    not_early <- outer(l2s$QTL_end+window_padding,l1s$QTL_start-window_padding,">=")
    not_late <- outer(l2s$QTL_start-window_padding,l1s$QTL_end+window_padding,"<=")
    overlap <- same_chr&not_late&not_early
    overlap <- mean(rowSums(overlap)>0)
    return(overlap)
  })
  return(c(mean(shared_a_b),mean(shared_b_a)))
}))

qtl_overlap <- data.frame(layer_comparisons[,1:2],qtl_overlap,qtl_overlap_perm)
colnames(qtl_overlap) <- c("A","B","A in B","B in A","A in B (permuted)","B in A (permuted)")
qtl_overlap[qtl_overlap=="pQTL"] <- "pQTL (sQTL)"
qtl_overlap[qtl_overlap=="pQTL2"] <- "pQTL (phosphoQTL)"
qtl_overlap[,3:6] <- round(qtl_overlap[,3:6],digits=2)
printRedmine(qtl_overlap)

genes <- intersect(eQTL$ensembl_gene_id,pQTL$ensembl_gene_id)


qtl_shared <- sapply(genes,FUN=function(g){
  e <- eQTL[eQTL$ensembl_gene_id==g,,drop=F]
  p <- pQTL[pQTL$ensembl_gene_id==g,,drop=F]
  same_chr <- outer(e$QTL_chr,p$QTL_chr,"==")
  not_early <- outer(e$QTL_end,p$QTL_start,">=")
  not_late <- outer(e$QTL_start,p$QTL_end,"<=")
  sum(same_chr&not_early&not_late)
})

qtl_shared_perm <- sapply(genes,FUN=function(g){
  perms <- sapply(1:10,FUN=function(i){
    g2 <- sample(genes,size = 1)
    e <- eQTL[eQTL$ensembl_gene_id==g,,drop=F]
    p <- pQTL[pQTL$ensembl_gene_id==g2,,drop=F]
    same_chr <- outer(e$QTL_chr,p$QTL_chr,"==")
    not_early <- outer(e$QTL_end,p$QTL_start,">=")
    not_late <- outer(e$QTL_start,p$QTL_end,"<=")
    sum(same_chr&not_early&not_late)
  })
  return(perms)
})


#make overview tables with descriptive qtl statistics
qtl_types <- names(legacy_results)
n_traits_ext <- c(n_traits[-2],phrQTL=879)
names(n_traits_ext)[2] <- "pQTL"
overview_tables <- lapply(list(legacy_results,swiq_results),FUN=function(tlist){
  out <- sapply(qtl_types,FUN=function(qtl_type){
    qtl_table <- tlist[[qtl_type]]
    n_t <- n_traits_ext[qtl_type]
    n_QTL <- nrow(qtl_table)
    n_QTL <- n_QTL/n_t
    pep <- "peptide"%in%colnames(qtl_table)
    if(!pep){
      traits_QTL <- length(unique(qtl_table$ensembl_gene_id))
    }else{
      traits_QTL <- length(unique(qtl_table$peptide))
    }
    prop_no_qtl <- 1-(traits_QTL/n_t)
    if(pep){
      n_local <- length(unique(qtl_table[qtl_table$local,"peptide"]))
    }else{
      n_local <- length(unique(qtl_table[qtl_table$local,"ensembl_gene_id"]))
    }
    prop_local <- n_local/n_t
    return(c(n_QTL,
             prop_no_qtl,
             prop_local))
  })
  rownames(out) <- c("number of QTL per trait",
                     "proportion of traits without QTL",
                     "proportion of traits with local QTL")
  return(out)
})
overview_tables <- abind(overview_tables,along = 3)
dimnames(overview_tables)[[3]] <- c("legacy","SWIQ")

#n qtl
plist <- lapply(qtl_types,FUN=function(qtl_type){
  pl <- gg_barplot(overview_tables[1,qtl_type,],
                   bar_color = c("grey","goldenrod"),
                   ylab = "# of QTL per trait",
                   title = qtl_type,
                   xlab = "",
                   bar_labels = dimnames(overview_tables)[[3]],axis_text_size = 25,title_size = 30,ylab_size = 25)
})
pl <- grid.arrange(grobs=plist, ncol = 2)
ggsave(filename = paste0(plot_path,"n_qtl.png"),device = "png",plot = pl)

#proportion of traits without a qtl
plist <- lapply(qtl_types,FUN=function(qtl_type){
  pl <- gg_barplot(overview_tables[2,qtl_type,]*100,
                   bar_color = c("grey","goldenrod"),
                   ylab = "% of traits w/o QTL",
                   title = qtl_type,
                   xlab = "",
                   bar_labels = dimnames(overview_tables)[[3]],axis_text_size = 25,title_size = 30,ylab_size = 25)
})
pl <- grid.arrange(grobs=plist, ncol = 2)
ggsave(filename = paste0(plot_path,"traits_wo_qtl.png"),device = "png",plot = pl)

#proportion of traits with a local qtl
plist <- lapply(qtl_types,FUN=function(qtl_type){
  pl <- gg_barplot(overview_tables[3,qtl_type,]*100,
                   bar_color = c("grey","goldenrod"),
                   ylab = "% of traits w local QTL",
                   title = qtl_type,
                   xlab = "",
                   bar_labels = dimnames(overview_tables)[[3]],axis_text_size = 25,title_size = 30,ylab_size = 25)
})
pl <- grid.arrange(grobs=plist, ncol = 2)
ggsave(filename = paste0(plot_path,"traits_w_local_qtl.png"),device = "png",plot = pl)

#qtl per trait by method
pl <- gg_barplot(overview_tables[1,,1],
                 bar_color = c("grey"),
                 ylab = "# of QTL per trait",
                 title = "legacy method",
                 xlab = "",
                 bar_labels = qtl_types,axis_text_size = 25,title_size = 30,ylab_size = 25)
pl2 <- gg_barplot(overview_tables[1,,2],
                 bar_color = c("goldenrod"),
                 ylab = "# of QTL per trait",
                 title = "SWIQ",
                 xlab = "",
                 bar_labels = qtl_types,axis_text_size = 25,title_size = 30,ylab_size = 25)
plg <- grid.arrange(pl,pl2, ncol = 1)
ggsave(filename = paste0(plot_path,"n_qtl_2.png"),device = "png",plot = plg)

#number of expected qtl per trait
legacy_counts <- sapply(1:10,FUN=function(i){
  sum(table(legacy_results$eQTL$ensembl_gene_id)==i)
})
legacy_counts <- c(5429-length(unique(legacy_results$eQTL$ensembl_gene_id)),legacy_counts)
swiq_counts <- sapply(1:10,FUN=function(i){
  sum(table(eQTL$ensembl_gene_id)==i)
})
swiq_counts <- c(5429-length(unique(eQTL$ensembl_gene_id)),swiq_counts)
legacy_exp <- dbinom(x = c(0:10),size = nrow(legacy_results$eQTL),prob = 1/5429)*5429
swiq_exp <- dbinom(x = c(0:10),size = nrow(eQTL),prob = 1/5429)*5429
barplot(rbind(swiq_exp,swiq_counts),beside = T)
barplot(rbind(legacy_exp,legacy_counts),beside = T)
p1 <- gg_barplot_matrix(rbind(swiq_exp,swiq_counts),within_pair_gap = 0.5,between_pair_gap = 0.75,bar_color = c("blue","goldenrod"),pair_labels = 0:10,ylab = "number of traits",xlab = "number of QTL",title = "SWIQ",legend_labels = c("expectation","SWIQ"),legend_title = "")
p2 <- gg_barplot_matrix(rbind(legacy_exp,legacy_counts),within_pair_gap = 0.5,between_pair_gap = 0.75,bar_color = c("blue","darkgrey"),pair_labels = 0:10,ylab = "number of traits",xlab = "number of QTL",title = "legacy method",legend_labels = c("expectation","legacy"),legend_title = "")
pl <- grid.arrange(grobs=list(p2,p1),ncol=1)
ggsave(plot = pl,filename = paste0(plot_path,"n_qtl_exp.png"))
