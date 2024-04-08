## This code reproduces the analysis of RNA sequencing data and visualises the overlap of proteomic and RNA data from the article
## "Similar, but not the same: multi-omics comparison of human valve interstitial cells and osteoblast osteogenic differentiation expanded with an estimation of data-dependent and data-independent PASEF"
## Code author: Polina Kuchur

# the analysis requires the following libraries to be installed and imported
library(ggplot2)
library(dplyr)
library(gplots)
library(plotly)
library(tools)
library(circlize)
library(tidyverse)
library(stringr)
library(reshape2)
library(rstatix)
library(ggpubr)
library(ggforce)
library(org.Hs.eg.db)
library(DESeq2)
library(RColorBrewer)
library(VennDiagram)
library(ComplexHeatmap)
library(ReactomePA)
library(clusterProfiler)
library(edgeR)
library(readxl)

####---Part-1--- Transcriptome data analysis -----------####

# move to the working directory
setwd("./transcriptome/")

# upload file with RNA-seq counts
count_matrix <- as.matrix(read.csv("./vic_ost_counts.csv", sep = ";", row.names = 1))
head(count_matrix, 2)

# filtering low-counted genes
counts <- subset(count_matrix, rowSums(count_matrix) >= 10)

# convert gene ensembl names to gene symbols
genes <- rownames(counts)
symbols <- mapIds(org.Hs.eg.db, keys = genes,
                  column = c('SYMBOL'),
                  keytype = 'ENSEMBL')
genes_symbol <- data.frame(genes, symbols)

counter=1
for(gene in genes_symbol$symbols){
  if (is.na(gene)) {
    genes_symbol[counter,2] <- genes_symbol[counter,1]
  }
  counter = counter+1
}

# change ensg to gene symbols where possible
genes_symbol$symbols <- make.names(genes_symbol$symbols, unique = TRUE)
rownames(counts) <- genes_symbol$symbols
nrow(counts)

# save counts after filtration
#write.csv(counts, file="rna_vic_ost_counts_cutoff10.csv")


###---Part-1.1--- Compare OB undif against OB dif ------###

# upload file with RNA-seq counts
count_matrix <- as.matrix(read.csv("./vic_ost_counts.csv", sep = ";", row.names = 1))
head(count_matrix)

# upload metadata
coldata_file <- read_xlsx("./vic_ost_metadata.xlsx", sheet = 3)

# transcriptome metadata
coldata <- data.frame(
  sample = coldata_file$sample,
  replicates = coldata_file$replicate,
  status = coldata_file$status,
  batch = coldata_file$batch,
  row.names = "sample")

# make sample characteristics as factors
coldata$replicate <- as.factor(coldata$replicate)
coldata$status <- as.factor(coldata$status)
coldata$batch <- as.factor(coldata$batch)

# select OBs
count_matrix <- count_matrix[,rownames(coldata)]

# check the order of samples in metadata and count_matrix
all(rownames(coldata) %in% colnames(count_matrix))

# prepare data for DESeq2
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~batch+status)

# merge technical replicates
dds <- collapseReplicates(dds, dds$replicates)

# filtration of low-counted genes
dds <- dds[rowSums(counts(dds)) >= 10,]

# set the reference condition
dds$status <- relevel(dds$status, ref = "10d.dif")

# median of ratios normalization
#dds <- estimateSizeFactors(dds)
#sizeFactors(dds)

# add gene symbols to gene properties
genes <- rownames(dds)
symbols <- mapIds(org.Hs.eg.db, keys = genes,
                  column = c('SYMBOL'),
                  keytype = 'ENSEMBL')
genes_symbol <- data.frame(genes, symbols)

counter=1
for(gene in genes_symbol$symbols){
  if (is.na(gene)) {
    genes_symbol[counter,2] <- genes_symbol[counter,1]
  }
  counter = counter+1
}
# change ensg to gene symbols where possible
rownames(dds) <- genes_symbol$symbols
nrow(dds)


## Start DESeq2
dds <- DESeq(dds)
resultsNames(dds)  # see all comparisons

# separate DESeq results
res_ob_undif_dif10 <- results(dds, contrast=c("status","48h.undif","10d.dif"))

# explore results
summary(res_ob_undif_dif10)

# export DESeq results in the file
#write.csv(as.data.frame(res_ob_undif_dif10[order(res_ob_undif_dif10$padj),] ), file="./rna_ost_undif_vs_dif_deseq.csv")


###---Part-1.2--- Compare VIC undif against VIC dif ------###

# upload file with RNA-seq counts
count_matrix <- as.matrix(read.csv("./vic_ost_counts.csv", sep = ";", row.names = 1))
head(count_matrix)

# upload metadata
coldata_file <- read_xlsx("./vic_ost_metadata.xlsx", sheet = 2)

# transcriptome metadata
coldata <- data.frame(
  sample = coldata_file$sample,
  replicates = coldata_file$replicate,
  status = coldata_file$status,
  batch = coldata_file$batch,
  row.names = "sample")

# make sample characteristics as factors
coldata$replicate <- as.factor(coldata$replicate)
coldata$status <- as.factor(coldata$status)
coldata$batch <- as.factor(coldata$batch)

# select VICs
count_matrix <- count_matrix[,rownames(coldata)]

# check the order of samples in metadata and count_matrix
all(rownames(coldata) %in% colnames(count_matrix))

# prepare data for DESeq2
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~batch+status)

# merge technical replicates
dds <- collapseReplicates(dds, dds$replicates)

# filtration of low-counted genes
dds <- dds[rowSums(counts(dds)) >= 10,]

# set the reference condition
dds$status <- relevel(dds$status, ref = "10d.dif")

# median of ratios normalization
#dds <- estimateSizeFactors(dds)
#sizeFactors(dds)

# add gene symbols to gene properties
genes <- rownames(dds)
symbols <- mapIds(org.Hs.eg.db, keys = genes,
                  column = c('SYMBOL'),
                  keytype = 'ENSEMBL')
genes_symbol <- data.frame(genes, symbols)

counter=1
for(gene in genes_symbol$symbols){
  if (is.na(gene)) {
    genes_symbol[counter,2] <- genes_symbol[counter,1]
  }
  counter = counter+1
}
# change ensg to gene symbols where possible
rownames(dds) <- genes_symbol$symbols
nrow(dds)


## Start DESeq2
dds <- DESeq(dds)
resultsNames(dds)  # see all comparisons

# separate DESeq results
res_vic_undif_dif10 <- results(dds, contrast=c("status","48h.undif","10d.dif"))

# explore results
summary(res_vic_undif_dif10)

# export DESeq results in the file
#write.csv(as.data.frame(res_vic_undif_dif10[order(res_vic_undif_dif10$padj),] ), file="./rna_vic_undif_vs_dif_deseq.csv")


###---Part-1.3--- Compare OB undif against VIC undif ------###

# upload RNA-seq counts
count_matrix <- as.matrix(read.csv("./vic_ost_counts.csv", sep = ";", row.names = 1))
head(count_matrix, 2)

# upload metadata
coldata_file <- read_xlsx("./vic_ost_metadata.xlsx", sheet = 4)

coldata <- data.frame(
  sample = coldata_file$sample,
  replicates = coldata_file$replicate,
  source = coldata_file$source,
  batch = coldata_file$batch,
  row.names = "sample")

# make samples characteristics as factors
coldata$replicate <- as.factor(coldata$replicate)
coldata$source <- as.factor(coldata$source)
coldata$batch <- as.factor(coldata$batch)

# check the order of samples in metadata and count_matrix
all(rownames(coldata) %in% colnames(count_matrix))

# get subset of control (48h.undif) samples from count_matrix (all samples)
count_matrix <- count_matrix[,rownames(coldata)]

# prepare data for DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~source)

dds <- collapseReplicates(dds, dds$replicates)

# set the reference condition
dds$source <- relevel(dds$source, ref = "VIC")

# median of ratios normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# filtration of low-counted genes
dds <- dds[rowSums(counts(dds)) >= 10,]

# add gene symbols to gene properties
genes <- rownames(dds)
symbols <- mapIds(org.Hs.eg.db, keys = genes,
                  column = c('SYMBOL'),
                  keytype = 'ENSEMBL')
genes_symbol <- data.frame(genes, symbols)

counter=1
for(gene in genes_symbol$symbols){
  if (is.na(gene)) {
    genes_symbol[counter,2] <- genes_symbol[counter,1]
  }
  counter = counter+1
}

# change ensg to gene symbols where possible
rownames(dds) <- genes_symbol$symbols
nrow(dds)
 
## Start DESeq2 analysis
dds <- DESeq(dds)
resultsNames(dds)  # see all comparisons 

# extract DESeq2 results 
res <- results(dds, contrast=c("source","OB","VIC"))

# summary of differentially expressed genes
summary(res, alpha = 0.05)

# export DESeq2 results in the file
#write.csv(as.data.frame(res[order(abs(res$log2FoldChange), decreasing = TRUE),] ), file="./rna_ost_vs_vic_undif_deseq.csv")


###---Part-1.4--- Compare OB 10 days of dif against VIC 10 days of dif -----###

# upload RNA-seq counts
count_matrix <- as.matrix(read.csv("./vic_ost_counts.csv", sep = ";", row.names = 1))
head(count_matrix, 2)

# upload metadata
coldata_file <- read_xlsx("./vic_ost_metadata.xlsx", sheet = 5)

# transcriptome metadata
coldata <- data.frame(
  sample = coldata_file$sample,
  replicates = coldata_file$replicate,
  source = coldata_file$source,
  batch = coldata_file$batch,
  row.names = "sample")

# make samples characteristics as factors
coldata$replicate <- as.factor(coldata$replicate)
coldata$source <- as.factor(coldata$source)
coldata$batch <- as.factor(coldata$batch)

# check the order of samples in metadata and count_matrix
all(rownames(coldata) %in% colnames(count_matrix))

# get samples from common dataset based on their names
count_matrix <- count_matrix[,rownames(coldata)]

# prepare data for DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~source)

# merge replicates
dds <- collapseReplicates(dds, dds$replicates)

# set the reference condition
dds$source <- relevel(dds$source, ref = "VIC")

# median of ratios normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# filtration of low-counted genes
dds <- dds[rowSums(counts(dds)) >= 10,]

# add gene symbols to gene properties
genes <- rownames(dds)
symbols <- mapIds(org.Hs.eg.db, keys = genes,
                  column = c('SYMBOL'),
                  keytype = 'ENSEMBL')
genes_symbol <- data.frame(genes, symbols)

counter=1
for(gene in genes_symbol$symbols){
  if (is.na(gene)) {
    genes_symbol[counter,2] <- genes_symbol[counter,1]
  }
  counter = counter+1
}

# change ensg to gene symbols where possible
rownames(dds) <- genes_symbol$symbols
nrow(dds)

## Start DESeq2 analysis
dds <- DESeq(dds)
resultsNames(dds)  # see all comparisons 

# extract DESeq2 results 
res <- results(dds, contrast=c("source","OB","VIC"))

# summary of dif. expressed genes
summary(res)

# export DESeq2 results in the file
#write.csv(as.data.frame(res[order(abs(res$log2FoldChange), decreasing = TRUE),] ), file="./rna_ost_vs_vic_dif10d_deseq.csv")



###---Part-2--- Overlap of transcriptome and proteome data -----###

## Figure 2a - Overlap of DEGs (RNA-seq) and DEPs (DIA, DDA, DIA-ML)

# proteome dia
dia_csv <- read.csv("../data_prep_to_github/dia_imp.csv")
dia_gene <- dia_csv$X
# proteome dda
dda_csv <- read.csv("../data_prep_to_github/dda_imp.csv")
dda_gene <- dda_csv$X
# proteome dia-ml
diaml_csv <- read.csv("../data_prep_to_github/diaml_imp.csv")  
diaml_gene <- diaml_csv$X
# transcriptome
rna_csv <- read.csv("./rna_vic_ost_counts_cutoff10.csv")
rna_gene <- rna_csv$X

# visualize overlaps by Venn diagram
myCol <- brewer.pal(4, "Pastel2")
venn.diagram(
  x = list(dia_gene,
           dda_gene,
           diaml_gene,
           rna_gene),
  category.names = c("DIA", 
                     "DDA",
                     "DIA-ML",
                     "RNA-seq"),
  filename = './venn_dia_dda_diaml_rna.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1500 , 
  width = 1500 , 
  resolution = 500,
  compression = "lzw",
  
  # Circles
  lwd = 1,
  lty = 'solid',
  fill = myCol,
  
  # Numbers
  cex = .5,
  fontface = "bold",
  fontfamily = "sans",
  main.col = "black",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  print.mode = c("raw","percent"),
  sigdigs = 2,
  cat.fontfamily = "sans"
  )


## Figure 3b - Enrichment analysis of DEGs and DEPs upregulated in OB
# enrichment is performed for undifferentiated (control) and differentiated cells separately

## DEPs control
dia_contr <- read.csv("../DIA/dia_ost_contVSvic_cont.csv")
dia_contr_upreg <- subset(dia_contr, dia_contr$logFC > 1 & dia_contr$adj.P.Val < 0.05)
dia_contr_upreg <- dia_contr_upreg$X

dda_contr <- read.csv("../DDA/dda_ost_contVSvic_cont.csv")
dda_contr_upreg <- subset(dda_contr, dda_contr$logFC > 1 & dda_contr$adj.P.Val < 0.05)
dda_contr_upreg <- dda_contr_upreg$X

diaml_contr <- read.csv("../DIA-ML/diaml_ost_contVSvic_cont.csv")
diaml_contr_upreg <- subset(diaml_contr, diaml_contr$logFC > 1 & diaml_contr$adj.P.Val < 0.05)
diaml_contr_upreg <- diaml_contr_upreg$X

# overlap of DEPs (Prot.rel contr.)
intersect_contr_upreg <- Reduce(intersect,list(dia_contr_upreg,
                                               dda_contr_upreg,
                                               diaml_contr_upreg))

# combine all DEPs (Prot.w contr.)
union_contr_upreg <- Reduce(union,list(dia_contr_upreg,
                                       dda_contr_upreg,
                                       diaml_contr_upreg))

## DEGs control
rna_contr_csv <- read.csv("./rna_ost_vs_vic_undif_deseq.csv")
rna_contr_significant <- subset(rna_contr_csv, rna_contr_csv$padj < 0.05 & rna_contr_csv$log2FoldChange > 1)

rna_diff_csv <- read.csv("./rna_ost_vs_vic_dif10d_deseq.csv")
rna_diff_significant <- subset(rna_diff_csv, rna_diff_csv$padj < 0.05 & rna_diff_csv$log2FoldChange > 1)


## DEPs dif10d 
dia_diff <- read.csv("../DIA/dia_ost_diffVSvic_diff.csv")
dia_diff_upreg <- subset(dia_diff, dia_diff$logFC > 1 & dia_diff$adj.P.Val < 0.05)
dia_diff_upreg <- dia_diff_upreg$X

dda_diff <- read.csv("../DDA/dda_ost_diffVSvic_diff.csv")
dda_diff_upreg <- subset(dda_diff, dda_diff$logFC > 1 & dda_diff$adj.P.Val < 0.05)
dda_diff_upreg <- dda_diff_upreg$X

diaml_diff <- read.csv("../DIA-ML/diaml_ost_diffVSvic_diff.csv")
diaml_diff_upreg <- subset(diaml_diff, diaml_diff$logFC > 1 & diaml_diff$adj.P.Val < 0.05)
diaml_diff_upreg <- diaml_diff_upreg$X

# overlap of DEPs (Prot.rel dif.)
intersect_diff_upreg <- Reduce(intersect,list(dia_diff_upreg,
                                               dda_diff_upreg,
                                               diaml_diff_upreg))
# combine DEPs (Prot.w dif.)
union_diff_upreg <- Reduce(union,list(dia_diff_upreg,
                                       dda_diff_upreg,
                                       diaml_diff_upreg))

prot_w_contr <- union_contr_upreg
prot_w_diff <- union_diff_upreg

prot_rel_contr <- intersect_contr_upreg
prot_rel_diff <- intersect_diff_upreg

rna_contr <- rna_contr_significant$X
rna_diff <- rna_diff_significant$X


## Enrichment analysis

# convert gene names to entrez
prot_w_contr.entrez <- clusterProfiler::bitr(prot_w_contr,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
prot_rel_contr.entrez <- clusterProfiler::bitr(prot_rel_contr,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
rna_contr.entrez <- clusterProfiler::bitr(rna_contr,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
rna_diff.entrez <- clusterProfiler::bitr(rna_diff,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
prot_w_diff.entrez <- clusterProfiler::bitr(prot_w_diff,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
prot_rel_diff.entrez <- clusterProfiler::bitr(prot_rel_diff,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

# create clusters
clusters <- list(prot_w_contr = prot_w_contr.entrez$ENTREZID,
                 prot_w_diff = prot_w_diff.entrez$ENTREZID,
                 prot_rel_contr = prot_rel_contr.entrez$ENTREZID,
                 prot_rel_diff = prot_rel_diff.entrez$ENTREZID,
                 rna_contr = rna_contr.entrez$ENTREZID,
                 rna_diff = rna_diff.entrez$ENTREZID)

## enrichGO
ck <- compareCluster(geneCluster = clusters, fun = enrichGO, OrgDb=org.Hs.eg.db, ont = "BP")

# make name of one position shorter
ck@compareClusterResult[["Description"]][[13]] <- "regulation of transmembrane receptor protein ser/thr kinase pathway"

# visualize these results
p1 <- dotplot(ck, 
              showCategory = 10,
              by = "GeneRatio",
              color = "p.adjust",
              font.size = 12.5,
              title="OB vs VIC upreg (Biological processes)",
              includeAll=TRUE) + scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=12))

# save the figure
# tiff(file="ostVSvic_biological_processes.tiff",
#      units = "in",
#      width = 7,
#      height = 15,
#      res=300)
gridExtra::grid.arrange(p1, ncol = 1)
dev.off()


## Figure 3c - Enrichment analysis of DEGs and DEPs upregulated in VIC
# enrichment is performed for undifferentiated (control) and differentiated cells separately

## DEPs control
dia_contr <- read.csv("../DIA/dia_ost_contVSvic_cont.csv")
dia_contr_upreg_vic <- subset(dia_contr, dia_contr$logFC < -1 & dia_contr$adj.P.Val < 0.05)
dia_contr_upreg_vic <- dia_contr_upreg_vic$X

dda_contr <- read.csv("../DDA/dda_ost_contVSvic_cont.csv")
dda_contr_upreg_vic <- subset(dda_contr, dda_contr$logFC < -1 & dda_contr$adj.P.Val < 0.05)
dda_contr_upreg_vic <- dda_contr_upreg_vic$X

diaml_contr <- read.csv("../DIA-ML/diaml_ost_contVSvic_cont.csv")
diaml_contr_upreg_vic <- subset(diaml_contr, diaml_contr$logFC < -1 & diaml_contr$adj.P.Val < 0.05)
diaml_contr_upreg_vic <- diaml_contr_upreg_vic$X

# intersection (Prot. rel contr.)
intersect_contr_upreg_vic <- Reduce(intersect,list(dia_contr_upreg_vic,
                                               dda_contr_upreg_vic,
                                               diaml_contr_upreg_vic))

union_contr_upreg_vic <- Reduce(union,list(dia_contr_upreg_vic,
                                       dda_contr_upreg_vic,
                                       diaml_contr_upreg_vic))
# DEGs
rna_contr_csv <- read.csv("./rna_ost_vs_vic_undif_deseq.csv")
rna_contr_significant_vic <- subset(rna_contr_csv, rna_contr_csv$padj < 0.05 & rna_contr_csv$log2FoldChange < -1)

rna_diff_csv <- read.csv("./rna_ost_vs_vic_dif10d_deseq.csv")
rna_diff_significant_vic <- subset(rna_diff_csv, rna_diff_csv$padj < 0.05 & rna_diff_csv$log2FoldChange < -1)

## DEPs diff 
dia_diff <- read.csv("../DIA/dia_ost_diffVSvic_diff.csv")
dia_diff_upreg_vic <- subset(dia_diff, dia_diff$logFC < -1 & dia_diff$adj.P.Val < 0.05)
dia_diff_upreg_vic <- dia_diff_upreg_vic$X

dda_diff <- read.csv("../DDA/dda_ost_diffVSvic_diff.csv")
dda_diff_upreg_vic <- subset(dda_diff, dda_diff$logFC < -1 & dda_diff$adj.P.Val < 0.05)
dda_diff_upreg_vic <- dda_diff_upreg_vic$X

diaml_diff <- read.csv("../DIA-ML/diaml_ost_diffVSvic_diff.csv")
diaml_diff_upreg_vic <- subset(diaml_diff, diaml_diff$logFC < -1 & diaml_diff$adj.P.Val < 0.05)
diaml_diff_upreg_vic <- diaml_diff_upreg_vic$X

# Prot.rel. diff.
intersect_diff_upreg_vic <- Reduce(intersect,list(dia_diff_upreg_vic,
                                              dda_diff_upreg_vic,
                                              diaml_diff_upreg_vic))

# Prot.w. diff.
union_diff_upreg_vic <- Reduce(union,list(dia_diff_upreg_vic,
                                      dda_diff_upreg_vic,
                                      diaml_diff_upreg_vic))

prot_w_contr <- union_contr_upreg_vic
prot_w_diff <- union_diff_upreg_vic

prot_rel_contr <- intersect_contr_upreg_vic
prot_rel_diff <- intersect_diff_upreg_vic

rna_contr <- rna_contr_significant_vic$X
rna_diff <- rna_diff_significant_vic$X


## Enrichment analysis

# convert gene names to entrez
prot_w_contr.entrez <- clusterProfiler::bitr(prot_w_contr,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
prot_rel_contr.entrez <- clusterProfiler::bitr(prot_rel_contr,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
rna_contr.entrez <- clusterProfiler::bitr(rna_contr,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
rna_diff.entrez <- clusterProfiler::bitr(rna_diff,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
prot_w_diff.entrez <- clusterProfiler::bitr(prot_w_diff,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
prot_rel_diff.entrez <- clusterProfiler::bitr(prot_rel_diff,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

# comparative clusters
clusters <- list(prot_w_contr = prot_w_contr.entrez$ENTREZID,
                 prot_w_diff = prot_w_diff.entrez$ENTREZID,
                 prot_rel_contr = prot_rel_contr.entrez$ENTREZID,
                 prot_rel_diff = prot_rel_diff.entrez$ENTREZID,
                 rna_contr = rna_contr.entrez$ENTREZID,
                 rna_diff = rna_diff.entrez$ENTREZID)

# enrichGO
ck <- compareCluster(geneCluster = clusters, fun = enrichGO, OrgDb=org.Hs.eg.db, ont = "BP")


p1 <- dotplot(ck, 
              showCategory = 10,
              by = "GeneRatio",
              color = "p.adjust",
              font.size = 12.5,
              title="VIC vs OB (Biological processes)",
              includeAll=TRUE) + scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=12))

# tiff(file="vicVSost_biological_processes.tiff.tiff",
#      units = "in",
#      width = 7,
#      height = 7,
#      res=300)
gridExtra::grid.arrange(p1, ncol = 1)
dev.off()



## Figure 2d - Heatmaps for Prot.rel (DIA, DDA & DIAML overlap) in control (48h.undif) and osteogenic differetiation (Ost. dif)

# control DEPs
dia_contr <- read.csv("../DIA/dia_ost_contVSvic_cont.csv")
dia_contr <- subset(dia_contr, abs(dia_contr$logFC) > 1.0 & dia_contr$adj.P.Val<0.05)
dia_contr <- dia_contr[,c(1,2)]

dda_contr <- read.csv("../DDA/dda_ost_contVSvic_cont.csv")
dda_contr <- subset(dda_contr, abs(dda_contr$logFC) > 1.0 & dda_contr$adj.P.Val<0.05)
dda_contr <- dda_contr[,c(1,2)]

diaml_contr <- read.csv("../DIA-ML/diaml_ost_contVSvic_cont.csv")
diaml_contr <- subset(diaml_contr, abs(diaml_contr$logFC) > 1.0 & diaml_contr$adj.P.Val<0.05)
diaml_contr <- diaml_contr[,c(1,2)]

df_list_control <- list(dia_contr, dda_contr, diaml_contr)

expr_control <- Reduce(function(x, y) merge(x, y, by = "X", all=FALSE), df_list_control)

colnames(expr_control) <- c("Gene", "dia_contr", "dda_contr", "diaml_contr")
rownames(expr_control) <- expr_control$Gene   

mat = as.matrix(expr_control[, grep("d", colnames(expr_control))])
mat_scaled = t(apply(mat, 1, scale))
base_mean = rowMeans(mat)
type = gsub("s\\d+_", "", colnames(mat))

# select colors
ha = HeatmapAnnotation(type = type, annotation_name_side = "left", col = list(type = c("dda_contr" = '#CBD5E8', "dia_contr" = '#B3E2CD', "diaml_contr" = '#FDCDAC')))
colors = gplots::colorpanel(256,"blue2","white","red2")

#pdf('./heatmap_ost_contr_vs_vic_control_dda_dia_diaml_overlap.pdf', width=5.5, height=8)
Heatmap(mat, name = "Expression", col = colors, row_title = NULL, show_row_dend = FALSE, show_row_names = FALSE, top_annotation = ha) + Heatmap(base_mean, name = "base mean", width = unit(15, "mm"))
dev.off()


## DEPs differentiation
dia_diff <- read.csv("../DIA/dia_ost_diffVSvic_diff.csv")
dia_diff <- subset(dia_diff, abs(dia_diff$logFC) > 1.0 & dia_diff$adj.P.Val<0.05)
dia_diff <- dia_diff[,c(1,2)]

dda_diff <- read.csv("../DDA/dda_ost_diffVSvic_diff.csv")
dda_diff <- subset(dda_diff, abs(dda_diff$logFC) > 1.0 & dda_diff$adj.P.Val<0.05)
dda_diff <- dda_diff[,c(1,2)]

diaml_diff <- read.csv("../DIA-ML/diaml_ost_diffVSvic_diff.csv")
diaml_diff <- subset(diaml_diff, abs(diaml_diff$logFC) > 1.0 & diaml_diff$adj.P.Val<0.05)
diaml_diff <- diaml_diff[,c(1,2)]

df_list_diff <- list(dia_diff, dda_diff, diaml_diff)

expr_diff <- Reduce(function(x, y) merge(x, y, by = "X", all=FALSE), df_list_diff)
rownames(expr_diff) <- expr_diff$X
colnames(expr_diff) <- c("Gene", "dia_dif", "dda_dif", "diaml_dif")

mat = as.matrix(expr_diff[, grep("d", colnames(expr_diff))])
mat_scaled = t(apply(mat, 1, scale))
base_mean = rowMeans(mat)
type = gsub("s\\d+_", "", colnames(mat))

# select colors
ha = HeatmapAnnotation(type = type, annotation_name_side = "left", col = list(type = c("dda_dif" = '#CBD5E8', "dia_dif" = '#B3E2CD', "diaml_dif" = '#FDCDAC')))
colors = gplots::colorpanel(256,"blue2","white","red2")

#pdf('./heatmap_ost_dif_vs_vic_dif_dda_dia_diaml_overlap.pdf', width=5, height=8)
Heatmap(mat, name = "Expression", col = colors, row_title = NULL, show_row_dend = FALSE, show_row_names = FALSE, top_annotation = ha) +
  Heatmap(base_mean, name = "base mean", width = unit(15, "mm"))
dev.off()

## Supplementary Materials 2. Correlation of Log2 Fold Changes for differentially expressed genes/proteins by RNA-seq and DIA proteomics.

##-- 2.1 in comparison of control osteoblasts and control VICs
rna_vc_c <- data.frame(read.csv("./rna_ost_vs_vic_undif_deseq.csv", header = TRUE))
rownames(rna_vc_c) <- rna_vc_c[,1]
colnames(rna_vc_c)[3] <- "logFC"

dia_vc_c <- data.frame(read.csv("../DIA/dia_ost_contVSvic_cont.csv", header = TRUE))
rownames(dia_vc_c) <- dia_vc_c[,1]

common <- intersect(rna_vc_c$X, dia_vc_c$X)  
FCrna <- rna_vc_c[rna_vc_c$X %in% common, c('X', 'logFC')]
FCrna <- FCrna[!duplicated(FCrna$X), ]
FCrna <- arrange(FCrna, FCrna$X)
colnames(FCrna) <- c('Gene_ID', 'logFC_rna')

FCdia <- dia_vc_c[dia_vc_c$X %in% common, c('X', 'logFC')]
FCdia <- FCdia[!duplicated(FCdia$X), ]
FCdia <- arrange(FCdia, FCdia$X)
colnames(FCdia) <- c('Gene_ID', 'logFC_dia')

Dia_rna <- data.frame(FCrna, FCdia)
head(Dia_rna)

#tiff('./rna_dia_vics_ost_contr.tiff', units="in", width=16, height=11, res=350, compression = 'lzw')
ggplot(Dia_rna, aes(x=logFC_rna , y=logFC_dia)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold changes between RNA and DIA proteomics ost_vs_vics_contr") + labs(x = "RNA logFC") + labs(y = "DIA logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_rna>1.5 & logFC_dia>1.5,as.character(Gene_ID),'')),hjust=0.5,vjust=-0.5, size = 6) + 
  geom_text(aes(label=ifelse(logFC_rna<(-1.5) & logFC_dia<(-1.5),as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6)
dev.off()

cor(x = Dia_rna$logFC_rna, y = Dia_rna$logFC_dia, method = "pearson")


##-- 2.2 in comparison of diff osteoblasts and diff VICs
rna_vc_d <- data.frame(read.csv("./rna_ost_vs_vic_dif10d_deseq.csv", header = TRUE))
rownames(rna_vc_d) <- rna_vc_d[,1]
colnames(rna_vc_d)[3] <- "logFC"

dia_vc_d <- data.frame(read.csv("../DIA/dia_ost_diffVSvic_diff.csv", header = TRUE)) 
rownames(dia_vc_d) <- dia_vc_d[,1]

common <- intersect(rna_vc_d$X, dia_vc_d$X)  
FCrna <- rna_vc_d[rna_vc_d$X %in% common, c('X', 'logFC')]
FCrna <- FCrna[!duplicated(FCrna$X), ]
FCrna <- arrange(FCrna, FCrna$X)
colnames(FCrna) <- c('Gene_ID', 'logFC_rna')

FCdia <- dia_vc_d[dia_vc_d$X %in% common, c('X', 'logFC')]
FCdia <- FCdia[!duplicated(FCdia$X), ]
FCdia <- arrange(FCdia, FCdia$X)
colnames(FCdia) <- c('Gene_ID', 'logFC_dia')

Dia_rna <- data.frame(FCrna, FCdia)
head(Dia_rna)

#tiff('./rna_dia_vics_ost_dif.tiff', units="in", width=16, height=11, res=350, compression = 'lzw')
ggplot(Dia_rna, aes(x=logFC_rna , y=logFC_dia)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold change between RNA and DIA proteomics ost_vic_diff") + labs(x = "RNA logFC") + labs(y = "DIA logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_rna>1.5 & logFC_dia>1.5,as.character(Gene_ID),'')),hjust=0.5,vjust=-0.5, size = 6) + 
  geom_text(aes(label=ifelse(logFC_rna<(-1.5) & logFC_dia<(-1.5),as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6)
dev.off()

cor(x = Dia_rna$logFC_rna, y = Dia_rna$logFC_dia, method = "pearson")


##-- 2.3 in comparison of diff osteoblasts and control osteoblasts
rna_ost <- data.frame(read.csv("./rna_ost_undif_vs_dif_deseq.csv", header = TRUE)) 
rownames(rna_ost) <- rna_ost[,1]
colnames(rna_ost)[3] <- "logFC"

dia_ost <- data.frame(read.csv("../DIA/dif_dia_ost_contVSost_diff.csv", header = TRUE))
rownames(dia_ost) <- dia_ost[,1]

common <- intersect(rna_ost$X, dia_ost$X)  
FCrna <- rna_ost[rna_ost$X %in% common, c('X', 'logFC')]
FCrna <- FCrna[!duplicated(FCrna$X), ]
FCrna <- arrange(FCrna, FCrna$X)
colnames(FCrna) <- c('Gene_ID', 'logFC_rna')

FCdia <- dia_ost[dia_ost$X %in% common, c('X', 'logFC')]
FCdia <- FCdia[!duplicated(FCdia$X), ]
FCdia <- arrange(FCdia, FCdia$X)
colnames(FCdia) <- c('Gene_ID', 'logFC_dia')

Dia_rna <- data.frame(FCrna, FCdia)
head(Dia_rna)

#tiff('./rna_dia_ost.tiff', units="in", width=16, height=11, res=350, compression = 'lzw')
ggplot(Dia_rna, aes(x=logFC_rna , y=logFC_dia)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold change between RNA and DIA proteomics ost") + labs(x = "RNA logFC") + labs(y = "DIA logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_rna>1.5 & logFC_dia>1.5,as.character(Gene_ID),'')),hjust=0.5,vjust=-0.5, size = 6) + 
  geom_text(aes(label=ifelse(logFC_rna<(-1.5) & logFC_dia<(-1.5),as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6)
dev.off()

cor(x = Dia_rna$logFC_rna, y = Dia_rna$logFC_dia, method = "pearson")


##-- 2.4 in comparison of diff VICs and control VICs
rna_vic <- data.frame(read.csv("./rna_vic_undif_vs_dif_deseq.csv", header = TRUE))
rownames(rna_vic) <- rna_vic[,1]
colnames(rna_vic)[3] <- "logFC"

dia_vic <- data.frame(read.csv("../DIA/dia_vic_contVSvic_diff.csv", header = TRUE))
rownames(dia_vic) <- dia_vic[,1]

common <- intersect(rna_vic$X, dia_vic$X)  
FCrna <- rna_vic[rna_vic$X %in% common, c('X', 'logFC')]
FCrna <- FCrna[!duplicated(FCrna$X), ]
FCrna <- arrange(FCrna, FCrna$X)
colnames(FCrna) <- c('Gene_ID', 'logFC_rna')

FCdia <- dia_vic[dia_vic$X %in% common, c('X', 'logFC')]
FCdia <- FCdia[!duplicated(FCdia$X), ]
FCdia <- arrange(FCdia, FCdia$X)
colnames(FCdia) <- c('Gene_ID', 'logFC_dia')

Dia_rna <- data.frame(FCrna, FCdia)
head(Dia_rna)

#tiff('./rna_dia_vic.tiff', units="in", width=16, height=11, res=350, compression = 'lzw')
ggplot(Dia_rna, aes(x=logFC_rna , y=logFC_dia)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold change between RNA and DIA proteomics vic") + labs(x = "RNA logFC") + labs(y = "DIA logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_rna>1.5 & logFC_dia>1.5,as.character(Gene_ID),'')),hjust=0.5,vjust=-0.5, size = 6) + 
  geom_text(aes(label=ifelse(logFC_rna<(-1.5) & logFC_dia<(-1.5),as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6)
dev.off()

cor(x = Dia_rna$logFC_rna, y = Dia_rna$logFC_dia, method = "pearson")



## Supplementary Materials 3. Clustering of osteoblasts and human valve interstitial cells (VICs) samples by principal component analysis (PCA).

count_matrix <- as.matrix(read.csv("./vic_ost_counts.csv", sep = ";", row.names = 1))
head(count_matrix, 2)

#transcriptome metadata
coldata_file <- read_xlsx("./vic_ost_metadata.xlsx", sheet = 1)

#transcriptome metadata
coldata <- data.frame(
  sample = coldata_file$sample,
  replicates = coldata_file$replicate,
  status = coldata_file$status,
  source = coldata_file$source,
  batch = coldata_file$batch,
  row.names = "sample")

#check the order of samples in metadata and count_matrix
all(rownames(coldata) %in% colnames(count_matrix))

##Start analysis
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~source+status)

dds <- collapseReplicates(dds, dds$replicates)

#set the reference condition
dds$status <- relevel(dds$status, ref = "48h.undif")
dds$source <- relevel(dds$source, ref = "OB")

#median of ratios normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

#filtration of low-counted genes
dds <- dds[rowSums(counts(dds)) >= 10,]

#add gene symbols to gene properties
genes <- rownames(dds)
symbols <- mapIds(org.Hs.eg.db, keys = genes,
                  column = c('SYMBOL'),
                  keytype = 'ENSEMBL')
genes_symbol <- data.frame(genes, symbols)

counter=1
for(gene in genes_symbol$symbols){
  if (is.na(gene)) {
    genes_symbol[counter,2] <- genes_symbol[counter,1]
  }
  counter = counter+1
}

#change ensg to gene symbols where possible
rownames(dds) <- genes_symbol$symbols
nrow(dds)

## DESeq
dds <- DESeq(dds)
resultsNames(dds)  # see all comparisons

rlt <- rlog(dds)  #rlog transformation


# tiff(file="./pca_rna.tiff",
#      units = "in",
#      width = 10,
#      height = 10,
#      res=300)
pcaData <- plotPCA(rlt, intgroup=c("source"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
names <- pcaData[,5]
pcaData[,5] <- str_replace_all(names, ".d20", "")
pcaData[,5] <- str_replace_all(pcaData[,5], ".d8", "")
pcaData[,5] <- str_replace_all(pcaData[,5], ".d9", "")
pcaData[,5] <- str_replace_all(pcaData[,5], ".d506", "")
pcaData[,5] <- str_replace_all(pcaData[,5], ".d508", "")
pcaData[,5] <- str_replace_all(pcaData[,5], ".d509", "")
colnames(pcaData) <-c("PC1","PC2","group","condition","replicate")

ggplot(pcaData, aes(PC1, PC2, shape=replicate, color=replicate)) +
  geom_point(size=3) +
  geom_mark_ellipse(aes(color = replicate)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  scale_colour_manual(values=c("#ff8F00", "#378EB8", "#2daf7A", "#999999"))

dev.off()



## Supplementary Materials 4. Venn diagrams representing the overlap of differentially expressed genes and proteins between control and osteogenic differentiation in osteoblasts and VICs

# legend
# a - RNA: ost_contrVSvic_control ⋂ ost_diffVSvic_diff
# e - RNA: vic_contrVSost_control ⋂ vic_diffVSost_diff
# f - prot.w: ost_contrVSvic_control ⋂ ost_diffVSvic_diff
# g - prot.w: vic_contrVSost_control ⋂ vic_diffVSost_diff
# h - prot.rel: ost_contrVSvic_control ⋂ ost_diffVSvic_diff
# i - prot.rel: vic_contrVSost_control ⋂ vic_diffVSost_diff

# a
rna_contr_csv <- read.csv("./rna_ost_vs_vic_undif_deseq.csv")
rna_contr_significant <- subset(rna_contr_csv, rna_contr_csv$padj < 0.05 & rna_contr_csv$log2FoldChange > 1)
ost_contrVSvic_control_rna <- rna_contr_significant$X

rna_diff_csv <- read.csv("./rna_ost_vs_vic_dif10d_deseq.csv")
rna_diff_significant <- subset(rna_diff_csv, rna_diff_csv$padj < 0.05 & rna_diff_csv$log2FoldChange > 1)
ost_diffVSvic_diff_rna <- rna_diff_significant$X

myCol <- brewer.pal(7, "Set1")
venn.diagram(
  x = list(ost_contrVSvic_control_rna,
           ost_diffVSvic_diff_rna),
  category.names = c("Control", 
                     "Differentiation"),
  filename = './figure14a.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 700 , 
  width = 700 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 1,
  lty = 1,
  fill = c("#377EB8", "#FF7F00"),
  main.cex = 10,
  
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.pos = c(200, 10),
  cat.dist = c(0.05,0.05),
  cat.fontfamily = "sans"
)


# b
rna_contr_csv <- read.csv("./rna_ost_vs_vic_undif_deseq.csv")
rna_contr_significant <- subset(rna_contr_csv, rna_contr_csv$padj < 0.05 & rna_contr_csv$log2FoldChange < -1)
vic_contrVSost_control_rna <- rna_contr_significant$X

rna_diff_csv <- read.csv("./rna_ost_vs_vic_dif10d_deseq.csv")
rna_diff_significant <- subset(rna_diff_csv, rna_diff_csv$padj < 0.05 & rna_diff_csv$log2FoldChange < -1)
vic_diffVSost_diff_rna <- rna_diff_significant$X

venn.diagram(
  x = list(vic_contrVSost_control_rna,
           vic_diffVSost_diff_rna),
  category.names = c("Control", 
                     "Differentiation"),
  filename = './figure14b.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 700 , 
  width = 700 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 1,
  lty = 1,
  fill = c("#999999", "#4DAF4A"),
  main.cex = 10,
  
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.pos = c(200, 20),
  cat.dist = c(0.05,0.05),
  cat.fontfamily = "sans"
)


# c
dia_contr <- read.csv("../DIA/dia_ost_contVSvic_cont.csv")
dia_contr_upreg <- subset(dia_contr, dia_contr$logFC > 1 & dia_contr$adj.P.Val < 0.05)
dia_contr_upreg <- dia_contr_upreg$X

dia_diff <- read.csv("../DIA/dia_ost_diffVSvic_diff.csv")
dia_diff_upreg <- subset(dia_diff, dia_diff$logFC > 1 & dia_diff$adj.P.Val < 0.05)
dia_diff_upreg <- dia_diff_upreg$X

dda_contr <- read.csv("../DDA/dda_ost_contVSvic_cont.csv")
dda_contr_upreg <- subset(dda_contr, dda_contr$logFC > 1 & dda_contr$adj.P.Val < 0.05)
dda_contr_upreg <- dda_contr_upreg$X

dda_diff <- read.csv("../DDA/dda_ost_diffVSvic_diff.csv")
dda_diff_upreg <- subset(dda_diff, dda_diff$logFC > 1 & dda_diff$adj.P.Val < 0.05)
dda_diff_upreg <- dda_diff_upreg$X

diaml_contr <- read.csv("../DIA-ML/diaml_ost_contVSvic_cont.csv")
diaml_contr_upreg <- subset(diaml_contr, diaml_contr$logFC > 1 & diaml_contr$adj.P.Val < 0.05)
diaml_contr_upreg <- diaml_contr_upreg$X

diaml_diff <- read.csv("../DIA-ML/diaml_ost_diffVSvic_diff.csv")
diaml_diff_upreg <- subset(diaml_diff, diaml_diff$logFC > 1 & diaml_diff$adj.P.Val < 0.05)
diaml_diff_upreg <- diaml_diff_upreg$X

#union of three proteomes
ost_contrVSvic_control_prot_w <- Reduce(union,list(dia_contr_upreg,
                                                   dda_contr_upreg,
                                                   diaml_contr_upreg))
ost_diffVSvic_diff_prot_w <- Reduce(union,list(dia_diff_upreg,
                                                       dda_diff_upreg,
                                                       diaml_diff_upreg))

venn.diagram(
  x = list(ost_contrVSvic_control_prot_w,
           ost_diffVSvic_diff_prot_w),
  category.names = c("Control", 
                     "Differentiation"),
  filename = './figure14c.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 700 , 
  width = 700 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 1,
  lty = 1,
  fill = c("#377EB8", "#FF7F00"),
  main.cex = 10,
  
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.pos = c(200, 20),
  cat.dist = c(0.05,0.05),
  cat.fontfamily = "sans"
)


# d
dia_contr <- read.csv("../DIA/dia_ost_contVSvic_cont.csv")
dia_contr_upreg_vic <- subset(dia_contr, dia_contr$logFC < -1 & dia_contr$adj.P.Val < 0.05)
dia_contr_upreg_vic <- dia_contr_upreg_vic$X

dia_diff <- read.csv("../DIA/dia_ost_diffVSvic_diff.csv")
dia_diff_upreg_vic <- subset(dia_diff, dia_diff$logFC < -1 & dia_diff$adj.P.Val < 0.05)
dia_diff_upreg_vic <- dia_diff_upreg_vic$X

dda_contr <- read.csv("../DDA/dda_ost_contVSvic_cont.csv")
dda_contr_upreg_vic <- subset(dda_contr, dda_contr$logFC < -1 & dda_contr$adj.P.Val < 0.05)
dda_contr_upreg_vic <- dda_contr_upreg_vic$X

dda_diff <- read.csv("../DDA/dda_ost_diffVSvic_diff.csv")
dda_diff_upreg_vic <- subset(dda_diff, dda_diff$logFC < -1 & dda_diff$adj.P.Val < 0.05)
dda_diff_upreg_vic <- dda_diff_upreg_vic$X

diaml_contr <- read.csv("../DIA-ML/diaml_ost_contVSvic_cont.csv")
diaml_contr_upreg_vic <- subset(diaml_contr, diaml_contr$logFC < -1 & diaml_contr$adj.P.Val < 0.05)
diaml_contr_upreg_vic <- diaml_contr_upreg_vic$X

diaml_diff <- read.csv("../DIA-ML/diaml_ost_diffVSvic_diff.csv")
diaml_diff_upreg_vic <- subset(diaml_diff, diaml_diff$logFC < -1 & diaml_diff$adj.P.Val < 0.05)
diaml_diff_upreg_vic <- diaml_diff_upreg_vic$X

#union of three proteomes
vic_contrVSost_control_prot_w <- Reduce(union,list(dia_contr_upreg_vic,
                                                   dda_contr_upreg_vic,
                                                   diaml_contr_upreg_vic))
vic_diffVSost_diff_prot_w <- Reduce(union,list(dia_diff_upreg_vic,
                                               dda_diff_upreg_vic,
                                               diaml_diff_upreg_vic))

venn.diagram(
  x = list(vic_contrVSost_control_prot_w,
           vic_diffVSost_diff_prot_w),
  category.names = c("Control", 
                     "Differentiation"),
  filename = './figure14d.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 700 , 
  width = 700 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 1,
  lty = 1,
  fill = c("#999999", "#4DAF4A"),
  main.cex = 10,
  
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.pos = c(200, 20),
  cat.dist = c(0.05,0.05),
  cat.fontfamily = "sans"
)


# e
dia_contr <- read.csv("../DIA/dia_ost_contVSvic_cont.csv")
dia_contr_upreg <- subset(dia_contr, dia_contr$logFC > 1 & dia_contr$adj.P.Val < 0.05)
dia_contr_upreg <- dia_contr_upreg$X

dia_diff <- read.csv("../DIA/dia_ost_diffVSvic_diff.csv")
dia_diff_upreg <- subset(dia_diff, dia_diff$logFC > 1 & dia_diff$adj.P.Val < 0.05)
dia_diff_upreg <- dia_diff_upreg$X

dda_contr <- read.csv("../DDA/dda_ost_contVSvic_cont.csv")
dda_contr_upreg <- subset(dda_contr, dda_contr$logFC > 1 & dda_contr$adj.P.Val < 0.05)
dda_contr_upreg <- dda_contr_upreg$X

dda_diff <- read.csv("../DDA/dda_ost_diffVSvic_diff.csv")
dda_diff_upreg <- subset(dda_diff, dda_diff$logFC > 1 & dda_diff$adj.P.Val < 0.05)
dda_diff_upreg <- dda_diff_upreg$X

diaml_contr <- read.csv("../DIA-ML/diaml_ost_contVSvic_cont.csv")
diaml_contr_upreg <- subset(diaml_contr, diaml_contr$logFC > 1 & diaml_contr$adj.P.Val < 0.05)
diaml_contr_upreg <- diaml_contr_upreg$X

diaml_diff <- read.csv("../DIA-ML/diaml_ost_diffVSvic_diff.csv")
diaml_diff_upreg <- subset(diaml_diff, diaml_diff$logFC > 1 & diaml_diff$adj.P.Val < 0.05)
diaml_diff_upreg <- diaml_diff_upreg$X

#intersect of three proteomes
ost_contrVSvic_control_prot_rel <- Reduce(intersect,list(dia_contr_upreg,
                                                   dda_contr_upreg,
                                                   diaml_contr_upreg))
ost_diffVSvic_diff_prot_rel <- Reduce(intersect,list(dia_diff_upreg,
                                               dda_diff_upreg,
                                               diaml_diff_upreg))

venn.diagram(
  x = list(ost_contrVSvic_control_prot_rel,
           ost_diffVSvic_diff_prot_rel),
  category.names = c("Control", 
                     "Differentiation"),
  filename = './figure14e.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 700 , 
  width = 700 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 1,
  lty = 1,
  fill = c("#377EB8", "#FF7F00"),
  main.cex = 10,
  
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.pos = c(200, 20),
  cat.dist = c(0.05,0.05),
  cat.fontfamily = "sans"
)


# f
dia_contr <- read.csv("../DIA/dia_ost_contVSvic_cont.csv")
dia_contr_upreg_vic <- subset(dia_contr, dia_contr$logFC < -1 & dia_contr$adj.P.Val < 0.05)
dia_contr_upreg_vic <- dia_contr_upreg_vic$X

dia_diff <- read.csv("../DIA/dia_ost_diffVSvic_diff.csv")
dia_diff_upreg_vic <- subset(dia_diff, dia_diff$logFC < -1 & dia_diff$adj.P.Val < 0.05)
dia_diff_upreg_vic <- dia_diff_upreg_vic$X

dda_contr <- read.csv("../DDA/dda_ost_contVSvic_cont.csv")
dda_contr_upreg_vic <- subset(dda_contr, dda_contr$logFC < -1 & dda_contr$adj.P.Val < 0.05)
dda_contr_upreg_vic <- dda_contr_upreg_vic$X

dda_diff <- read.csv("../DDA/dda_ost_diffVSvic_diff.csv")
dda_diff_upreg_vic <- subset(dda_diff, dda_diff$logFC < -1 & dda_diff$adj.P.Val < 0.05)
dda_diff_upreg_vic <- dda_diff_upreg_vic$X

diaml_contr <- read.csv("../DIA-ML/diaml_ost_contVSvic_cont.csv")
diaml_contr_upreg_vic <- subset(diaml_contr, diaml_contr$logFC < -1 & diaml_contr$adj.P.Val < 0.05)
diaml_contr_upreg_vic <- diaml_contr_upreg_vic$X

diaml_diff <- read.csv("../DIA-ML/diaml_ost_diffVSvic_diff.csv")
diaml_diff_upreg_vic <- subset(diaml_diff, diaml_diff$logFC < -1 & diaml_diff$adj.P.Val < 0.05)
diaml_diff_upreg_vic <- diaml_diff_upreg_vic$X

#intersect three proteomes
vic_contrVSost_control_prot_rel <- Reduce(intersect,list(dia_contr_upreg_vic,
                                                   dda_contr_upreg_vic,
                                                   diaml_contr_upreg_vic))
vic_diffVSost_diff_prot_rel <- Reduce(intersect,list(dia_diff_upreg_vic,
                                               dda_diff_upreg_vic,
                                               diaml_diff_upreg_vic))

venn.diagram(
  x = list(vic_contrVSost_control_prot_rel,
           vic_diffVSost_diff_prot_rel),
  category.names = c("Control", 
                     "Differentiation"),
  filename = './figure14f.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 700 , 
  width = 700 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 1,
  lty = 1,
  fill = c("#999999", "#4DAF4A"),
  main.cex = 10,
  
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.pos = c(200, 20),
  cat.dist = c(0.05,0.05),
  cat.fontfamily = "sans"
)


## Supplementary Materials 5. Boxplots for selective proteins are exclusively upregulated during VICs' osteogenic differentiation compared to osteoblasts differentiation.

dia <- read.csv("../data_prep_to_github/dia_toNA.csv", sep = ",")
gene = "GLIPR2"  # change to MAOA, AKAP2 or IRAG1
gene_subset <- subset(dia, dia$X == gene) 

# remove control samples
gene_subset <- gene_subset %>% select(-matches("contr"))

t_gene_subset <- t(gene_subset)
t_gene_subset <- t_gene_subset[-1]
t_gene_subset_num <- transform(t_gene_subset, as.numeric(t_gene_subset))
value <- t_gene_subset_num$X_data

headers <- str_split_fixed(colnames(gene_subset), "_", 3)
headers <- headers[,1:2]
header_joined <- paste(headers[,1], headers[,2], sep="_")
header_joined <- header_joined[-1]

status <- list(header_joined)
status <- rapply(status,function(x) ifelse(x=="Ost_contr","Ost. control",x), how = "replace")
status <- rapply(status,function(x) ifelse(x=="Ost_dif","Ost. differentiation",x), how = "replace")
status <- rapply(status,function(x) ifelse(x=="VIC_contr","VIC control",x), how = "replace")
status <- rapply(status,function(x) ifelse(x=="VIC_dif","VIC differentiation",x), how = "replace")

dataframe <- data.frame(source = status,
                        counts = value)
colnames(dataframe) <- c("Source", "counts")

data_mod <- melt(dataframe, id.vars='Source', measure.vars=c('counts')) 
data_mod <- transform(data_mod, value = as.numeric(value))
colnames(data_mod)[3] <- "Counts"

# create boxplot
# tiff(file="./glipr2_dia_dif.tiff",
#      units = "in",
#      width = 5,
#      height = 5,
#      res=300)
ggplot(data_mod, aes(x = Source, y = Counts)) +
  geom_boxplot(aes(fill = Source), show.legend = FALSE) +
  scale_fill_manual(values=c("#FF7F00", "#4DAF4A")) +
  stat_compare_means(method = "t.test") +
  theme_classic2() +
  labs(title = "GLIPR2 (DIA-PASEF)")  # change title in accordance with gene
dev.off()
