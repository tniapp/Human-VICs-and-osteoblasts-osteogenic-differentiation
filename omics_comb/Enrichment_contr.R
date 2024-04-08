### This part includes analysis of DDA dataset
## Opening the data and performing VSN normalization (normalization was disabled in NAguideR in previous step)
#setwd("your directory")


### 1 vic_vs_ob_contr
contr_dda <- data.frame(read.csv("dda_ost_contVSvic_cont.csv", header = TRUE)) 
rownames(contr_dda) <- contr_dda[,1]

contr_dia <- data.frame(read.csv("dia_ost_contVSvic_cont.csv", header = TRUE)) 
rownames(contr_dia) <- contr_dia[,1]

contr_diaml <- data.frame(read.csv("diaml_ost_contVSvic_cont.csv", header = TRUE)) 
rownames(contr_diaml) <- contr_diaml[,1]

contr_tr <- data.frame(read.csv("RNA_vic_undif_vs_ob_undif_deg.csv", header = TRUE)) 
rownames(contr_tr) <- contr_tr[,1]
contr_tr$log2FoldChange <- contr_tr$log2FoldChange*-1


common <- intersect(contr_dia$X, contr_tr$X)  
FCdia <- contr_dia[contr_dia$X %in% common, c('X', 'logFC')]
FCdia <- FCdia[!duplicated(FCdia$X), ]
FCdia <- arrange(FCdia, FCdia$X)
colnames(FCdia) <- c('Gene_ID', 'logFC_dia')

FCtr <- contr_tr[contr_tr$X %in% common, c('X', 'log2FoldChange')]
FCtr <- FCtr[!duplicated(FCtr$X), ]
FCtr <- arrange(FCtr, FCtr$X)
colnames(FCtr) <- c('Gene_ID', 'logFC_tr')

Tr_dia <- data.frame(FCdia, FCtr)
head(Tr_dia)

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
head(contr_dda)
up.genes_contr_dda <- contr_dda[contr_dda$logFC > 1 & contr_dda$adj.P.Val < 0.05, 1] 
dn.genes_contr_dda <- contr_dda[contr_dda$logFC < -1 & contr_dda$adj.P.Val < 0.05, 1]
bkgd.genes_contr_dda <- contr_dda[,1]

up.genes.entrez_contr_dda <- clusterProfiler::bitr(up.genes_contr_dda,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
dn.genes.entrez_contr_dda <- clusterProfiler::bitr(dn.genes_contr_dda,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez_contr_dda <- clusterProfiler::bitr(bkgd.genes_contr_dda,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

head(contr_dia)
up.genes_contr_dia <- contr_dia[contr_dia$logFC > 1 & contr_dia$adj.P.Val < 0.05, 1] 
dn.genes_contr_dia <- contr_dia[contr_dia$logFC < -1 & contr_dia$adj.P.Val < 0.05, 1]
bkgd.genes_contr_dia <- contr_dia[,1]

up.genes.entrez_contr_dia <- clusterProfiler::bitr(up.genes_contr_dia,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
dn.genes.entrez_contr_dia <- clusterProfiler::bitr(dn.genes_contr_dia,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez_contr_dia <- clusterProfiler::bitr(bkgd.genes_contr_dia,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


head(contr_diaml)
up.genes_contr_diaml <- contr_diaml[contr_diaml$logFC > 1 & contr_diaml$adj.P.Val < 0.05, 1] 
dn.genes_contr_diaml <- contr_diaml[contr_diaml$logFC < -1 & contr_diaml$adj.P.Val < 0.05, 1]
bkgd.genes_contr_diaml <- contr_diaml[,1]

up.genes.entrez_contr_diaml <- clusterProfiler::bitr(up.genes_contr_diaml,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
dn.genes.entrez_contr_diaml <- clusterProfiler::bitr(dn.genes_contr_diaml,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez_contr_diaml <- clusterProfiler::bitr(bkgd.genes_contr_diaml,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


head(contr_tr)
up.genes_contr_tr <- contr_tr[contr_tr$log2FoldChange > 1 & contr_tr$padj < 0.05, 1] 
up.genes_contr_tr <- up.genes_contr_tr[1:86]
dn.genes_contr_tr <- contr_tr[contr_tr$log2FoldChange < -1 & contr_tr$padj < 0.05, 1]
dn.genes_contr_tr <- dn.genes_contr_tr[1:110]
bkgd.genes_contr_tr <- contr_tr[,1]

up.genes.entrez_contr_tr <- clusterProfiler::bitr(up.genes_contr_tr,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
dn.genes.entrez_contr_tr <- clusterProfiler::bitr(dn.genes_contr_tr,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez_contr_tr <- clusterProfiler::bitr(bkgd.genes_contr_tr,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

#upregulated

library(clusterProfiler )

egobp_dia <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez_contr_dia[[2]],
  universe = bkgd.genes.entrez_contr_dia[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

egobp_dda <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez_contr_dda[[2]],
  universe = bkgd.genes.entrez_contr_dda[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

egobp_diaml <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez_contr_diaml[[2]],
  universe = bkgd.genes.entrez_contr_diaml[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

egobp_tr <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez_contr_tr[[2]],
  universe = bkgd.genes.entrez_contr_tr[[2]],
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

library(dplyr)
egobp <- full_join(egobp_dia, egobp_diaml, by = "Description")
colnames(egobp) <- c("GO BP", "DIA", "DIA-ML")



library(ggplot2)

egobp[is.na(egobp)] <- 0

rownames(egobp) <- egobp[,1]
egobp <- egobp[,-1]

tiff('go bp.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
egobp %>%
  rownames_to_column(var = "GO_BP") %>%
  gather(Dataset, count, -GO_BP) %>%
  ggplot() +
  geom_point(aes(x = Dataset, y = GO_BP, size = count) )
dev.off()





#CC
egocc_dia <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez_contr_dia[[2]],
  universe = bkgd.genes.entrez_contr_dia[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

egocc_dda <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez_contr_dda[[2]],
  universe = bkgd.genes.entrez_contr_dda[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

egocc_diaml <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez_contr_diaml[[2]],
  universe = bkgd.genes.entrez_contr_diaml[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

egocc_tr <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez_contr_tr[[2]],
  universe = bkgd.genes.entrez_contr_tr[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)


egocc_dia <- data.frame(egocc_dia)
egocc_dda <- data.frame(egocc_dda)
egocc_diaml <- data.frame(egocc_diaml)
egocc_tr <- data.frame(egocc_tr)


egocc_dia <- egocc_dia[,c(2,9)]
egocc_diaml <- egocc_diaml[,c(2,9)]
egocc_tr <- egocc_tr[,c(2,9)]

egocc <- full_join(egocc_dia, egocc_diaml, by = "Description")
egocc <- full_join(egocc, egocc_tr, by = "Description")

colnames(egocc) <- c("GO CC", "DIA", "DIA-ML", "RNA-seq")


egocc[is.na(egocc)] <- 0

rownames(egocc) <- egocc[,1]
egocc <- egocc[,-1]

tiff('go cc.tiff', units="in", width=15, height=9, res=300, compression = 'lzw')
egocc %>%
  rownames_to_column(var = "GO_CC") %>%
  gather(Dataset, count, -GO_CC) %>%
  ggplot() +
  geom_point(aes(x = Dataset, y = GO_CC, size = count) )
dev.off()



#kegg
kegg_dia <- clusterProfiler::enrichKEGG(
  gene     = up.genes.entrez_contr_dia[[2]],
  universe = bkgd.genes.entrez_contr_dia[[2]],
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05)

kegg_dda <- clusterProfiler::enrichKEGG(
  gene     = up.genes.entrez_contr_dda[[2]],
  universe = bkgd.genes.entrez_contr_dda[[2]],
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05)

kegg_diaml <- clusterProfiler::enrichKEGG(
  gene     = up.genes.entrez_contr_diaml[[2]],
  universe = bkgd.genes.entrez_contr_diaml[[2]],
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05)

kegg_tr <- clusterProfiler::enrichKEGG(
  gene     = up.genes.entrez_contr_tr[[2]],
  universe = bkgd.genes.entrez_contr_tr[[2]],
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
  gene     = dn.genes.entrez_contr_dia[[2]],
  universe = bkgd.genes.entrez_contr_dia[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

degobp_dda <- clusterProfiler::enrichGO(
  gene     = dn.genes.entrez_contr_dda[[2]],
  universe = bkgd.genes.entrez_contr_dda[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

degobp_diaml <- clusterProfiler::enrichGO(
  gene     = dn.genes.entrez_contr_diaml[[2]],
  universe = bkgd.genes.entrez_contr_diaml[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

degobp_tr <- clusterProfiler::enrichGO(
  gene     = dn.genes.entrez_contr_tr[[2]],
  universe = bkgd.genes.entrez_contr_tr[[2]],
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

tiff('contr_down_go bp.tiff', units="in", width=15, height=14, res=300, compression = 'lzw')
degobp %>%
  rownames_to_column(var = "GO_BP") %>%
  gather(Dataset, count, -GO_BP) %>%
  ggplot() +
  geom_point(aes(x = Dataset, y = GO_BP, size = count) )
dev.off()


