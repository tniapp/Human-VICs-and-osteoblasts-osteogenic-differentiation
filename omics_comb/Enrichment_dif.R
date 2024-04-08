### This part includes analysis of DDA dataset
## Opening the data and performing VSN normalization (normalization was disabled in NAguideR in previous step)
#setwd("your directory")


### 1 vic_vs_ob_dif
dif_dda <- data.frame(read.csv("dda_ost_diffVSvic_diff.csv", header = TRUE)) 
rownames(dif_dda) <- dif_dda[,1]

dif_dia <- data.frame(read.csv("dia_ost_diffVSvic_diff.csv", header = TRUE)) 
rownames(dif_dia) <- dif_dia[,1]

dif_diaml <- data.frame(read.csv("diaml_ost_diffVSvic_diff.csv", header = TRUE)) 
rownames(dif_diaml) <- dif_diaml[,1]

dif_tr <- data.frame(read.csv("RNA_vic_dif_vs_ob_dif_deg.csv", header = TRUE)) 
rownames(dif_tr) <- dif_tr[,1]
dif_tr$log2FoldChange <- dif_tr$log2FoldChange*-1

library(dplyr)
common <- intersect(dif_dia$X, dif_tr$X)  
FCdia <- dif_dia[dif_dia$X %in% common, c('X', 'logFC')]
FCdia <- FCdia[!duplicated(FCdia$X), ]
FCdia <- arrange(FCdia, FCdia$X)
colnames(FCdia) <- c('Gene_ID', 'logFC_dia')

FCtr <- dif_tr[dif_tr$X %in% common, c('X', 'log2FoldChange')]
FCtr <- FCtr[!duplicated(FCtr$X), ]
FCtr <- arrange(FCtr, FCtr$X)
colnames(FCtr) <- c('Gene_ID', 'logFC_tr')

Tr_dia <- data.frame(FCdia, FCtr)
head(Tr_dia)

library(ggplot2)
ggplot(Tr_dia, aes(x=logFC_dia , y=logFC_tr)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold change between DIA and tr") + labs(x = "DIA logFC") + labs(y = "Tr logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_dia>2 & logFC_tr>2,as.character(Gene_ID),'')),hjust=0,vjust=-0.5) + 
  geom_text(aes(label=ifelse(logFC_dia<(-2) & logFC_tr<(-2),as.character(Gene_ID),'')),hjust=0,vjust=-0.5)




#Pathway enrichment analysis
library(rWikiPathways)

load.libs <- c(
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3")
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
library("org.Hs.eg.db")


#extract genes
head(dif_dda)
up.genes_dif_dda <- dif_dda[dif_dda$logFC > 1 & dif_dda$adj.P.Val < 0.05, 1] 
dn.genes_dif_dda <- dif_dda[dif_dda$logFC < -1 & dif_dda$adj.P.Val < 0.05, 1]
bkgd.genes_dif_dda <- dif_dda[,1]

up.genes.entrez_dif_dda <- clusterProfiler::bitr(up.genes_dif_dda,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
dn.genes.entrez_dif_dda <- clusterProfiler::bitr(dn.genes_dif_dda,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez_dif_dda <- clusterProfiler::bitr(bkgd.genes_dif_dda,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

head(dif_dia)
up.genes_dif_dia <- dif_dia[dif_dia$logFC > 1 & dif_dia$adj.P.Val < 0.05, 1] 
dn.genes_dif_dia <- dif_dia[dif_dia$logFC < -1 & dif_dia$adj.P.Val < 0.05, 1]
bkgd.genes_dif_dia <- dif_dia[,1]

up.genes.entrez_dif_dia <- clusterProfiler::bitr(up.genes_dif_dia,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
dn.genes.entrez_dif_dia <- clusterProfiler::bitr(dn.genes_dif_dia,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez_dif_dia <- clusterProfiler::bitr(bkgd.genes_dif_dia,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


head(dif_diaml)
up.genes_dif_diaml <- dif_diaml[dif_diaml$logFC > 1 & dif_diaml$adj.P.Val < 0.05, 1] 
dn.genes_dif_diaml <- dif_diaml[dif_diaml$logFC < -1 & dif_diaml$adj.P.Val < 0.05, 1]
bkgd.genes_dif_diaml <- dif_diaml[,1]

up.genes.entrez_dif_diaml <- clusterProfiler::bitr(up.genes_dif_diaml,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
dn.genes.entrez_dif_diaml <- clusterProfiler::bitr(dn.genes_dif_diaml,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez_dif_diaml <- clusterProfiler::bitr(bkgd.genes_dif_diaml,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


head(dif_tr)
up.genes_dif_tr <- dif_tr[dif_tr$log2FoldChange > 1 & dif_tr$padj < 0.05, 1] 
up.genes_dif_tr <- up.genes_dif_tr[1:107]
dn.genes_dif_tr <- dif_tr[dif_tr$log2FoldChange < -1 & dif_tr$padj < 0.05, 1]
dn.genes_dif_tr <- dn.genes_dif_tr[1:124]
bkgd.genes_dif_tr <- dif_tr[,1]

up.genes.entrez_dif_tr <- clusterProfiler::bitr(up.genes_dif_tr,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
dn.genes.entrez_dif_tr <- clusterProfiler::bitr(dn.genes_dif_tr,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez_dif_tr <- clusterProfiler::bitr(bkgd.genes_dif_tr,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

#upregulated

library(clusterProfiler )

egobp_dia <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez_dif_dia[[2]],
  universe = bkgd.genes.entrez_dif_dia[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

egobp_dda <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez_dif_dda[[2]],
  universe = bkgd.genes.entrez_dif_dda[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

egobp_diaml <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez_dif_diaml[[2]],
  universe = bkgd.genes.entrez_dif_diaml[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

egobp_tr <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez_dif_tr[[2]],
  universe = bkgd.genes.entrez_dif_tr[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)


egobp_dia <- data.frame(egobp_dia)
egobp_dda <- data.frame(egobp_dda)
egobp_diaml <- data.frame(egobp_diaml)
egobp_tr <- data.frame(egobp_tr)


egobp_dia <- egobp_dia[,c(2,9)]
egobp_diaml <- egobp_diaml[,c(2,9)]
egobp_dda <- egobp_dda[,c(2,9)]
egobp_tr <- egobp_tr[,c(2,9)]

library(dplyr)
egobp <- full_join(egobp_dia, egobp_diaml, by = "Description")
egobp <- full_join(egobp, egobp_dda, by = "Description")
egobp <- full_join(egobp, egobp_tr, by = "Description")


colnames(egobp) <- c("GO BP", "DIA", "DIA-ML", "DDA", "RNA-seq")

egobp[is.na(egobp)] <- 0

rownames(egobp) <- egobp[,1]
egobp <- egobp[,-1]
library(tidyverse)
tiff('diff_go bp.tiff', units="in", width=15, height=20, res=300, compression = 'lzw')
egobp %>%
  rownames_to_column(var = "GO_BP") %>%
  gather(Dataset, count, -GO_BP) %>%
  ggplot() +
  geom_point(aes(x = Dataset, y = GO_BP, size = count) )
dev.off()


#kegg
kegg_dia <- clusterProfiler::enrichKEGG(
  gene     = up.genes.entrez_dif_dia[[2]],
  universe = bkgd.genes.entrez_dif_dia[[2]],
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05)

kegg_dda <- clusterProfiler::enrichKEGG(
  gene     = up.genes.entrez_dif_dda[[2]],
  universe = bkgd.genes.entrez_dif_dda[[2]],
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05)

dotplot(kegg_dda)
kegg_diaml <- clusterProfiler::enrichKEGG(
  gene     = up.genes.entrez_dif_diaml[[2]],
  universe = bkgd.genes.entrez_dif_diaml[[2]],
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05)

kegg_tr <- clusterProfiler::enrichKEGG(
  gene     = up.genes.entrez_dif_tr[[2]],
  universe = bkgd.genes.entrez_dif_tr[[2]],
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05)



kegg_dia <- data.frame(kegg_dia)
kegg_dda <- data.frame(kegg_dda)
kegg_diaml <- data.frame(kegg_diaml)
kegg_tr <- data.frame(kegg_tr)


kegg_dia <- kegg_dia[,c(2,9)]
kegg_diaml <- kegg_diaml[,c(2,9)]
kegg_dda <- kegg_dda[,c(2,9)]
kegg_tr <- kegg_tr[,c(2,9)]

KEGG <- full_join(kegg_dia, kegg_diaml, by = "Description")
KEGG <- full_join(KEGG, kegg_dda, by = "Description")
KEGG <- full_join(KEGG, kegg_tr, by = "Description")

colnames(KEGG) <- c("KEGG", "DIA", "DIA-ML", "DDA", "RNA-seq")


KEGG[is.na(KEGG)] <- 0

rownames(KEGG) <- KEGG[,1]
KEGG <- KEGG[,-1]

tiff('KEGG.tiff', units="in", width=8, height=6, res=300, compression = 'lzw')
KEGG %>%
  rownames_to_column(var = "KEGG") %>%
  gather(Dataset, count, -KEGG) %>%
  ggplot() +
  geom_point(aes(x = Dataset, y = KEGG, size = count) )
dev.off()



#downregulated

degobp_dia <- clusterProfiler::enrichGO(
  gene     = dn.genes.entrez_dif_dia[[2]],
  universe = bkgd.genes.entrez_dif_dia[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

degobp_dda <- clusterProfiler::enrichGO(
  gene     = dn.genes.entrez_dif_dda[[2]],
  universe = bkgd.genes.entrez_dif_dda[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

degobp_diaml <- clusterProfiler::enrichGO(
  gene     = dn.genes.entrez_dif_diaml[[2]],
  universe = bkgd.genes.entrez_dif_diaml[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "MF",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

degobp_tr <- clusterProfiler::enrichGO(
  gene     = dn.genes.entrez_dif_tr[[2]],
  universe = bkgd.genes.entrez_dif_tr[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)


degobp_dia <- data.frame(degobp_dia)
degobp_dda <- data.frame(degobp_dda)
degobp_diaml <- data.frame(degobp_diaml)
degobp_tr <- data.frame(degobp_tr)


degobp_dia <- degobp_dia[,c(2,9)]
degobp_diaml <- degobp_diaml[,c(2,9)]
degobp_dda <- degobp_dda[,c(2,9)]

library(dplyr)
degobp <- full_join(degobp_dia, degobp_diaml, by = "Description")
degobp <- full_join(degobp, degobp_dda, by = "Description")

colnames(degobp) <- c("GO BP", "DIA", "DIA-ML", "DDA")


degobp[is.na(degobp)] <- 0

rownames(degobp) <- degobp[,1]
degobp <- degobp[,-1]

tiff('dif_down_go bp.tiff', units="in", width=15, height=14, res=300, compression = 'lzw')
degobp %>%
  rownames_to_column(var = "GO_BP") %>%
  gather(Dataset, count, -GO_BP) %>%
  ggplot() +
  geom_point(aes(x = Dataset, y = GO_BP, size = count) )
dev.off()


