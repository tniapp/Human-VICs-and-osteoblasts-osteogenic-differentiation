# opening the data - using gene-level interference to reduce level of noise between different softwares and to easy compare our data with RNA-seq

# unzip report..pr_matrix_DiA_ph archive before run this code
setwd("./qualitative_analysis/")

library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(RColorBrewer)
library(VennDiagram)

###-----------------------------------------------------###
###-------------------- Control ------------------------###
###-----------------------------------------------------###

pept_dia <- data.frame(read.table("../data_prep_to_github/report..pr_matrix_DiA_ph/report..pr_matrix_DiA_ph.tsv", sep = '\t', header = TRUE))
prot_dia <- data.frame(read.table("../data_prep_to_github/report..gg_matrix_DiA_ph.tsv", sep = '\t', header = TRUE))

unique_pept_dia <- unique(pept_dia[,c('Genes','Stripped.Sequence')])
dia_un <- data.frame(table(unique_pept_dia$Genes))
dia_names <- dia_un[dia_un$Freq>1,]
dia_names <- as.character(dia_names$Var1)
prot_dia_f <- prot_dia[prot_dia$Genes %in% dia_names,]

head(prot_dia_f)

rownames(prot_dia_f) = make.names(prot_dia_f[,1], unique=TRUE)
prot_dia_f <- prot_dia_f[,-1]


OB_contr <- prot_dia_f[,c(1,3,5,7,9,11)]
head(OB_contr)

VIC_contr <- prot_dia_f[,c(13,15,17,19,21,23)]
head(VIC_contr)

OB_contr1 <- OB_contr[rowSums(is.na(OB_contr)) < 2,]
VIC_contr1 <- VIC_contr[rowSums(is.na(VIC_contr)) < 2,]

OB_contr_genes <- rownames(OB_contr1)
VIC_contr_genes <- rownames(VIC_contr1)

OB_contr_genes <-  sub("\\..*$", "", OB_contr_genes)
VIC_contr_genes <-  sub("\\..*$", "", VIC_contr_genes)

venn.diagram(
  x = list(OB_contr_genes,
           VIC_contr_genes),
  category.names = c("OB DIA", 
                     "VIC DIA"),
  filename = './venn_dia_ob_vic_all_control.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 1,
  lty = 1,
  fill = c("#0073C2FF", "grey"),
  main.cex = 10,
  
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.pos = c(220, 140),
  cat.dist = c(0.05,0.05),
  cat.fontfamily = "sans"
)

#write.csv(VIC_contr_genes, "VIC_contr_genes_DIA.csv")
#write.csv(OB_contr_genes, "OB_contr_genes_DIA.csv")

# select genes which are unique for OB or for VIC
OB_contr_genes_unique <- setdiff(OB_contr_genes, VIC_contr_genes)
#write.csv(OB_contr_genes_unique, "OB_contr_genes_unique_DIA.csv")

VIC_contr_genes_unique <- setdiff(VIC_contr_genes, OB_contr_genes)
#write.csv(VIC_contr_genes_unique, "VIC_contr_genes_unique_DIA.csv")

##-------------------- OB DIA --------------------------
## enrichment for unique genes of OB

OB_contr_genes_unique.entrez <- clusterProfiler::bitr(OB_contr_genes_unique,
                                                  fromType = "SYMBOL",
                                                  toType = "ENTREZID", 
                                                  OrgDb = org.Hs.eg.db)
OB_contr_genes_unique.entrez <- unique(OB_contr_genes_unique.entrez)

## enrichGO
OB_dia_unique_control_go <- enrichGO(OB_contr_genes_unique.entrez$ENTREZID,
                                     "org.Hs.eg.db",
                                     keyType = "ENTREZID",
                                     ont = "BP",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 0.2,
                                     minGSSize = 10,
                                     maxGSSize = 500)

# tiff(file="./enrichGO_ob_dia_unique_genes.tiff",
#           units = "in",
#           width = 6,
#           height = 6,
#           res=300)
dotplot(OB_dia_unique_control_go,
        showCategory = 10,
        color = "p.adjust",
        font.size = 15,
        title="OB unique DIA | GO: Biological processes") +
  scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=13), plot.title = element_text(hjust = 0.5))
dev.off()


##--------------------- VIC DIA --------------------------
## enrichment for unique genes of VIC

VIC_contr_genes_unique.entrez <- clusterProfiler::bitr(VIC_contr_genes_unique,
                                                      fromType = "SYMBOL",
                                                      toType = "ENTREZID", 
                                                      OrgDb = org.Hs.eg.db)
VIC_contr_genes_unique.entrez <- unique(VIC_contr_genes_unique.entrez)

## enrichGO
VIC_dia_unique_control_go <- enrichGO(VIC_contr_genes_unique.entrez$ENTREZID,
                                     "org.Hs.eg.db",
                                     keyType = "ENTREZID",
                                     ont = "BP",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 0.2,
                                     minGSSize = 10,
                                     maxGSSize = 500)

# tiff(file="./enrichGO_vic_dia_unique_genes.tiff",
#           units = "in",
#           width = 6,
#           height = 5,
#           res=300,
#      compression = "lzw")
dotplot(VIC_dia_unique_control_go,
        showCategory = 20,
        color = "p.adjust",
        font.size = 15,
        title="VIC unique DIA | GO: Biological processes") +
  scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=13), plot.title = element_text(hjust = 0.5))
dev.off()


##----------------------- OB RNA --------------------------

rna_filtered <- as.matrix(read.csv("../transcriptome/vic_ost_counts.csv", sep = ";", row.names = 1))
rna_filtered[rna_filtered == 0] = NA

rna_OB_contr <- rna_filtered[,c(1,2,4,5,8,9)]
head(rna_OB_contr)

rna_VIC_contr <- rna_filtered[,c(11,12,16,17,19,20)]
head(rna_VIC_contr)

rna_OB_contr1 <- rna_OB_contr[rowSums(is.na(rna_OB_contr)) < 2 & rowSums(rna_OB_contr, na.rm = TRUE) >= 10,]
rna_VIC_contr1 <- rna_VIC_contr[rowSums(is.na(rna_VIC_contr)) < 2 & rowSums(rna_VIC_contr, na.rm = TRUE) >= 10,]

rna_OB_contr_genes <- rownames(rna_OB_contr1) ## ENSG format
rna_VIC_contr_genes <- rownames(rna_VIC_contr1)

#write.csv(rna_OB_contr_genes, "rna_OB_contr_genes.csv")
#write.csv(rna_VIC_contr_genes, "rna_VIC_contr_genes.csv")

venn.diagram(
  x = list(rna_OB_contr_genes,
           rna_VIC_contr_genes),
  category.names = c("OB RNA", 
                     "VIC RNA"),
  filename = './venn_rna_ob_vic_all_control.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 1,
  lty = 1,
  fill = c("#0073C2FF", "grey"),
  main.cex = 10,
  
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.pos = c(220, 140),
  cat.dist = c(0.05,0.05),
  cat.fontfamily = "sans"
)

# select unique genes
rna_OB_contr_genes_unique <- setdiff(rna_OB_contr_genes,rna_VIC_contr_genes)
rna_VIC_contr_genes_unique <- setdiff(rna_VIC_contr_genes,rna_OB_contr_genes)


##----------------------- OB RNA  --------------------------
## enrichment for unique RNA genes of OB

rna_OB_contr_genes_unique.entrez <- clusterProfiler::bitr(rna_OB_contr_genes_unique,
                                                      fromType = "ENSEMBL",
                                                      toType = "ENTREZID", 
                                                      OrgDb = org.Hs.eg.db)
rna_OB_contr_genes_unique.entrez <- unique(rna_OB_contr_genes_unique.entrez)

## enrichGO
OB_rna_unique_control_go <- enrichGO(rna_OB_contr_genes_unique.entrez$ENTREZID,
                                     "org.Hs.eg.db",
                                     keyType = "ENTREZID",
                                     ont = "BP",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 0.2,
                                     minGSSize = 10,
                                     maxGSSize = 500)

# tiff(file="./enrichGO_ob_rna_unique_genes.tiff",
#           units = "in",
#           width = 6,
#           height = 6,
#           res=300, compression = "lzw")
dotplot(OB_rna_unique_control_go,
        showCategory = 10,
        color = "p.adjust",
        font.size = 15,
        title="OB unique RNA-seq | GO: Biological processes") +
  scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=13), plot.title = element_text(hjust = 0.5))
dev.off()


##----------------------- VIC RNA  --------------------------
## enrichment for unique genes RNA-seq of VIC

rna_VIC_contr_genes_unique.entrez <- clusterProfiler::bitr(rna_VIC_contr_genes_unique,
                                                       fromType = "ENSEMBL",
                                                       toType = "ENTREZID", 
                                                       OrgDb = org.Hs.eg.db)
rna_VIC_contr_genes_unique.entrez <- unique(rna_VIC_contr_genes_unique.entrez)

## enrichGO
VIC_rna_unique_control_go <- enrichGO(rna_VIC_contr_genes_unique.entrez$ENTREZID,
                                      "org.Hs.eg.db",
                                      keyType = "ENTREZID",
                                      ont = "BP",
                                      pvalueCutoff = 0.05,
                                      pAdjustMethod = "BH",
                                      qvalueCutoff = 0.2,
                                      minGSSize = 10,
                                      maxGSSize = 500)

# tiff(file="./enrichGO_vic_rna_unique_genes.tiff",
#      units = "in",
#      width = 6,
#      height = 6.5,
#      res=300, 
#      compression = "lzw")
dotplot(VIC_rna_unique_control_go,
        showCategory = 10,
        color = "p.adjust",
        font.size = 15,
        title="VIC unique RNA-seq | GO: Biological processes") +
  scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=13), plot.title = element_text(hjust = 0.5))
dev.off()


###-----------------------------------------------------###
###----------------- Differentiation -------------------###
###-----------------------------------------------------###

##------------------- Proteome -------------------

pept_dia <- data.frame(read.table("../data_prep_to_github/report..pr_matrix_DiA_ph/report..pr_matrix_DiA_ph.tsv", sep = '\t', header = TRUE))
prot_dia <- data.frame(read.table("../data_prep_to_github/report..gg_matrix_DiA_ph.tsv", sep = '\t', header = TRUE))

unique_pept_dia <- unique(pept_dia[,c('Genes','Stripped.Sequence')])
dia_un <- data.frame(table(unique_pept_dia$Genes))
dia_names <- dia_un[dia_un$Freq>1,]
dia_names <- as.character(dia_names$Var1)
prot_dia_f <- prot_dia[prot_dia$Genes %in% dia_names,]

head(prot_dia_f)

rownames(prot_dia_f) = make.names(prot_dia_f[,1], unique=TRUE)
prot_dia_f <- prot_dia_f[,-1]


OB_diff <- prot_dia_f[,c(2,4,6,8,10,12)]
head(OB_diff)

VIC_diff <- prot_dia_f[,c(14,16,18,20,22,24)]
head(VIC_diff)

OB_diff1 <- OB_diff[rowSums(is.na(OB_diff)) < 2,]
VIC_diff1 <- VIC_diff[rowSums(is.na(VIC_diff)) < 2,]

OB_diff_genes <- rownames(OB_diff1)
VIC_diff_genes <- rownames(VIC_diff1)

OB_diff_genes <-  sub("\\..*$", "", OB_diff_genes)
VIC_diff_genes <-  sub("\\..*$", "", VIC_diff_genes)


venn.diagram(
  x = list(OB_diff_genes,
           VIC_diff_genes),
  category.names = c("OB dif DIA", 
                     "VIC dif DIA"),
  filename = './venn_dia_ob_vic_all_diff.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 1,
  lty = 1,
  fill = c("#0073C2FF", "grey"),
  main.cex = 10,
  
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.pos = c(220, 140),
  cat.dist = c(0.04,0.04),
  cat.fontfamily = "sans"
)

# select genes which are unique for OB or for VIC
OB_dif_genes_unique <- setdiff(OB_diff_genes, VIC_diff_genes)
#write.csv(OB_dif_genes_unique, "OB_dif_genes_unique_DIA.csv")

VIC_dif_genes_unique <- setdiff(VIC_diff_genes, OB_diff_genes)
#write.csv(VIC_dif_genes_unique, "VIC_dif_genes_unique_DIA.csv")


##---------------------- OB DIA ------------------------
## enrichment for unique genes of OB
OB_dif_genes_unique.entrez <- clusterProfiler::bitr(OB_dif_genes_unique,
                                                      fromType = "SYMBOL",
                                                      toType = "ENTREZID", 
                                                      OrgDb = org.Hs.eg.db)
OB_dif_genes_unique.entrez <- unique(OB_dif_genes_unique.entrez)

## enrichGO
OB_dia_unique_dif_go <- enrichGO(OB_dif_genes_unique.entrez$ENTREZID,
                                     "org.Hs.eg.db",
                                     keyType = "ENTREZID",
                                     ont = "BP",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 0.2,
                                     minGSSize = 10,
                                     maxGSSize = 500)

# tiff(file="./enrichGO_ob_dia_unique_genes_diff.tiff",
#           units = "in",
#           width = 6,
#           height = 7,
#           res=300)
dotplot(OB_dia_unique_dif_go,
        showCategory = 10,
        color = "p.adjust",
        font.size = 15,
        title="OB unique DIA dif | GO: Biological processes") +
  scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=13), plot.title = element_text(hjust = 0.5))
dev.off()


##--------------------- VIC DIA -----------------------
## enrichment for unique genes of VIC

VIC_dif_genes_unique.entrez <- clusterProfiler::bitr(VIC_dif_genes_unique,
                                                       fromType = "SYMBOL",
                                                       toType = "ENTREZID", 
                                                       OrgDb = org.Hs.eg.db)
VIC_dif_genes_unique.entrez <- unique(VIC_dif_genes_unique.entrez)

## enrichGO
VIC_dia_unique_dif_go <- enrichGO(VIC_dif_genes_unique.entrez$ENTREZID,
                                      "org.Hs.eg.db",
                                      keyType = "ENTREZID",
                                      ont = "BP",
                                      pvalueCutoff = 0.05,
                                      pAdjustMethod = "BH",
                                      qvalueCutoff = 0.2,
                                      minGSSize = 10,
                                      maxGSSize = 500)

# tiff(file="./enrichGO_vic_dia_unique_genes_dif.tiff",
#           units = "in",
#           width = 6,
#           height = 5,
#           res=300,
#      compression = "lzw")
dotplot(VIC_dia_unique_dif_go,
        showCategory = 20,
        color = "p.adjust",
        font.size = 15,
        title="VIC unique DIA dif | GO: Biological processes") +
  scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=13), plot.title = element_text(hjust = 0.5))
dev.off()


##------------------- Transcriptome -------------------

rna_filtered <- as.matrix(read.csv("../transcriptome/vic_ost_counts.csv", sep = ";", row.names = 1))
rna_filtered[rna_filtered == 0] = NA

rna_OB_dif <- rna_filtered[,c(3,6,7,10)]
head(rna_OB_dif)

rna_VIC_dif <- rna_filtered[,c(13,14,15,18,21,22)]
head(rna_VIC_dif)

rna_OB_dif1 <- rna_OB_dif[rowSums(is.na(rna_OB_dif)) < 2 & rowSums(rna_OB_dif, na.rm = TRUE) >= 10,]
rna_VIC_dif1 <- rna_VIC_dif[rowSums(is.na(rna_VIC_dif)) < 2 & rowSums(rna_VIC_dif, na.rm = TRUE) >= 10,]

rna_OB_dif_genes <- rownames(rna_OB_dif1) ## ENSG format
rna_VIC_dif_genes <- rownames(rna_VIC_dif1)

#write.csv(rna_OB_dif_genes, "rna_OB_dif_genes.csv")
#write.csv(rna_VIC_dif_genes, "rna_VIC_dif_genes.csv")

venn.diagram(
  x = list(rna_OB_dif_genes,
           rna_VIC_dif_genes),
  category.names = c("OB dif RNA", 
                     "VIC dif RNA"),
  filename = './venn_rna_ob_vic_all_dif.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 1,
  lty = 1,
  fill = c("#0073C2FF", "grey"),
  main.cex = 10,
  
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.pos = c(200, 145),
  cat.dist = c(0.05,0.05),
  cat.fontfamily = "sans"
)

# select unique genes
rna_OB_dif_genes_unique <- setdiff(rna_OB_dif_genes,rna_VIC_dif_genes)
rna_VIC_dif_genes_unique <- setdiff(rna_VIC_dif_genes,rna_OB_dif_genes)


##----------------------- OB RNA  --------------------------
## enrichment for unique RNA genes of OB

rna_OB_dif_genes_unique.entrez <- clusterProfiler::bitr(rna_OB_dif_genes_unique,
                                                          fromType = "ENSEMBL",
                                                          toType = "ENTREZID", 
                                                          OrgDb = org.Hs.eg.db)
rna_OB_dif_genes_unique.entrez <- unique(rna_OB_dif_genes_unique.entrez)

## enrichGO
OB_rna_unique_dif_go <- enrichGO(rna_OB_dif_genes_unique.entrez$ENTREZID,
                                     "org.Hs.eg.db",
                                     keyType = "ENTREZID",
                                     ont = "BP",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 0.2,
                                     minGSSize = 10,
                                     maxGSSize = 500)

# tiff(file="./enrichGO_ob_rna_unique_genes_dif.tiff",
#           units = "in",
#           width = 6,
#           height = 6,
#           res=300, compression = "lzw")
dotplot(OB_rna_unique_dif_go,
        showCategory = 10,
        color = "p.adjust",
        font.size = 15,
        title="OB unique RNA-seq dif | GO: Biological processes") +
  scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=13), plot.title = element_text(hjust = 0.5))
dev.off()


##----------------------- VIC RNA  --------------------------
## enrichment for unique genes RNA-seq of VIC

rna_VIC_dif_genes_unique.entrez <- clusterProfiler::bitr(rna_VIC_dif_genes_unique,
                                                           fromType = "ENSEMBL",
                                                           toType = "ENTREZID", 
                                                           OrgDb = org.Hs.eg.db)
rna_VIC_dif_genes_unique.entrez <- unique(rna_VIC_dif_genes_unique.entrez)

## enrichGO
VIC_rna_unique_dif_go <- enrichGO(rna_VIC_dif_genes_unique.entrez$ENTREZID,
                                      "org.Hs.eg.db",
                                      keyType = "ENTREZID",
                                      ont = "BP",
                                      pvalueCutoff = 0.05,
                                      pAdjustMethod = "BH",
                                      qvalueCutoff = 0.2,
                                      minGSSize = 10,
                                      maxGSSize = 500)

# tiff(file="./enrichGO_vic_rna_unique_genes_dif.tiff",
#      units = "in",
#      width = 6,
#      height = 6.5,
#      res=300,
#      compression = "lzw")
dotplot(VIC_rna_unique_dif_go,
        showCategory = 10,
        color = "p.adjust",
        font.size = 15,
        title="VIC unique RNA-seq dif | GO: Biological processes") +
  scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=13), plot.title = element_text(hjust = 0.5))
dev.off()
