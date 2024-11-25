### This code includes data preparation and qualititative data comparison steps. Here we used outputs of specific software or pipelines to obtain proteomic or transcriptomic raw data.
### We have four datasets: ddaPASEF proteomics (FragPipe), diaPASEF proteomics with empirical spectral library (FragPip + DIA-NN), diaPASEF with ML-generated spectral library (DIA-NN), transcriptomic data (STAR)
#setwd("your directory")

## overlap in identified peptides between the three protromics datasets
# opening the data
pept_dia <- data.frame(read.table("report..pr_matrix_DiA_ph.tsv", sep = '\t', header = TRUE))
pept_diaml <- data.frame(read.table("report..pr_matrix_diaml.tsv", sep = '\t', header = TRUE))
pept_dda <- data.frame(read.table("combined_peptide_DDA.tsv", sep = '\t', header = TRUE))

library(VennDiagram)
library(RColorBrewer)
myCol2 <- brewer.pal(3, "Pastel2")

# Venn diagrams
venn.diagram(
  x = list(pept_dia$Stripped.Sequence, pept_diaml$Stripped.Sequence, pept_dda$Peptide.Sequence),
  category.names = c("DIA" , "DIA-ML" , "DDA"),
  filename = '#venn_diagramm_pept.png',
  resolution = 300,
  fill = myCol2,
  print.mode=c("raw","percent"),
  output=F
)

venn.diagram(
  x = list(pept_dia$Stripped.Sequence, pept_dda$Peptide.Sequence),
  category.names = c("DIA" , "DDA"),
  filename = '#venn_diagramm_pept_dda_dia.png',
  resolution = 300,
  print.mode=c("raw","percent"),
  output=F
)

venn.diagram(
  x = list(pept_dia$Stripped.Sequence,  pept_diaml$Stripped.Sequence),
  category.names = c("DIA" , "DIA-ML"),
  filename = '#venn_diagramm_pept_dia_diaml.png',
  resolution = 300,
  print.mode=c("raw","percent"),
  output=F
)

venn.diagram(
  x = list(pept_dia$Genes, pept_diaml$Genes, pept_dda$Gene),
  category.names = c("DIA" , "DIA-ML" , "DDA"),
  filename = '#venn_diagramm_pept_genes.png',
  resolution = 300,
  fill = myCol2,
  print.mode=c("raw","percent"),
  output=F
)


## preparation of the three protromics datasets for further analysis and demosntration of thier overlap on protein level
# opening the data - using gene-level interference to reduce level of noise between different softwares and to easy compare our data with RNA-seq
prot_diaml <- data.frame(read.table("report..gg_matrix_diaml.tsv", sep = '\t', header = TRUE))
prot_dia <- data.frame(read.table("report..gg_matrix_DiA_ph.tsv", sep = '\t', header = TRUE))
library(readxl)
prot_dda <- data.frame(read_excel("prot_dda.xlsx", sheet = 1))

#removing contaminants
prot_dda_f <- prot_dda[!grepl("contam", prot_dda$Protein), ]
prot_dda_cont <- prot_dda[grepl("contam", prot_dda$Protein), ]

rownames(prot_dia) = make.names(prot_dia[,1], unique=TRUE)
prot_dia <- prot_dia[,-1]
prot_dia_f <- prot_dia[!(row.names(prot_dia) %in% prot_dda_cont$Gene), ]

rownames(prot_diaml) <- prot_diaml[,1]
prot_diaml <- prot_diaml[,-1]
prot_diaml_f <- prot_diaml[!(row.names(prot_diaml) %in% prot_dda_cont$Gene), ]

#remplace NA in gene name in DDA data
library(dplyr)
prot_dda_f$Gene <- ifelse(is.na(prot_dda_f$Gene), prot_dda_f$Entry.Name, prot_dda_f$Gene)


venn.diagram(
  x = list(rownames(prot_dia_f), rownames(prot_diaml_f), prot_dda_f$Gene),
  category.names = c("DIA" , "DIA-ML" , "DDA"),
  filename = '#venn_diagramm_prot_raw.png',
  resolution = 300,
  fill = myCol2,
  print.mode=c("raw","percent"),
  output=F
)


#removing proteins with less than 2 peptides
prot_dda_f <- prot_dda_f[!prot_dda_f$Combined.Total.Peptides == 1, ]
tail(prot_dda_f)
rownames(prot_dda_f) = make.names(prot_dda_f[,4], unique=TRUE)
prot_dda_f <- prot_dda_f[,111:134]

unique_pept_dia <- unique(pept_dia[,c('Genes','Stripped.Sequence')])
dia_un <- data.frame(table(unique_pept_dia$Genes))
dia_names <- dia_un[dia_un$Freq>1,]
dia_names <- as.character(dia_names$Var1)
prot_dia_f <- prot_dia[dia_names,]


unique_pept_diaml <- unique(pept_diaml[,c('Genes','Stripped.Sequence')])
ml_un <- data.frame(table(unique_pept_diaml$Genes))
ml_names <- ml_un[ml_un$Freq>1,]
ml_names <- as.character(ml_names$Var1)
prot_diaml_f <- prot_diaml[ml_names,]

#changing sample names
library(readxl)
rename_dda <- data.frame(read_excel("rename.xlsx", sheet = 1))
colnames(prot_dda_f) <- rename_dda[,2]

rename_dia <- data.frame(read_excel("rename.xlsx", sheet = 2))
colnames(prot_dia_f) <- rename_dia[,2]
colnames(prot_diaml_f) <- rename_dia[,2]

#Removing VICs donor number 426
prot_dda_f <- prot_dda_f[,-c(14, 20)]
prot_dia_f <- prot_dia_f[,-c(15, 16)]
prot_diaml_f <- prot_diaml_f[,-c(15, 16)]

str(prot_dda_f)
prot_dda_f$Ost_contr_11 <- as.numeric(prot_dda_f$Ost_contr_11)
prot_dda_f$Ost_contr_12 <- as.numeric(prot_dda_f$Ost_contr_12)
prot_dda_f$Ost_contr_20 <- as.numeric(prot_dda_f$Ost_contr_20)
prot_dda_f$Ost_contr_22 <- as.numeric(prot_dda_f$Ost_contr_22)
prot_dda_f$Ost_contr_8 <- as.numeric(prot_dda_f$Ost_contr_8)
prot_dda_f$Ost_contr_9 <- as.numeric(prot_dda_f$Ost_contr_9)
prot_dda_f$Ost_dif_11 <- as.numeric(prot_dda_f$Ost_dif_11)
prot_dda_f$Ost_dif_12 <- as.numeric(prot_dda_f$Ost_dif_12)
prot_dda_f$Ost_dif_20 <- as.numeric(prot_dda_f$Ost_dif_20)
prot_dda_f$Ost_dif_22 <- as.numeric(prot_dda_f$Ost_dif_22)
prot_dda_f$Ost_dif_8 <- as.numeric(prot_dda_f$Ost_dif_8)
prot_dda_f$Ost_dif_9 <- as.numeric(prot_dda_f$Ost_dif_9)
prot_dda_f$VIC_contr_213 <- as.numeric(prot_dda_f$VIC_contr_213)
prot_dda_f$VIC_contr_431 <- as.numeric(prot_dda_f$VIC_contr_431)
prot_dda_f$VIC_contr_432 <- as.numeric(prot_dda_f$VIC_contr_432)
prot_dda_f$VIC_contr_482 <- as.numeric(prot_dda_f$VIC_contr_482)
prot_dda_f$VIC_contr_494 <- as.numeric(prot_dda_f$VIC_contr_494)
prot_dda_f$VIC_dif_213 <- as.numeric(prot_dda_f$VIC_dif_213)
prot_dda_f$VIC_dif_431 <- as.numeric(prot_dda_f$VIC_dif_431)
prot_dda_f$VIC_dif_432 <- as.numeric(prot_dda_f$VIC_dif_432)
prot_dda_f$VIC_dif_482 <- as.numeric(prot_dda_f$VIC_dif_482)
prot_dda_f$VIC_dif_494 <- as.numeric(prot_dda_f$VIC_dif_494)

boxplot(prot_dda_f, outline = TRUE, col = cols, main = "Raw data", names=colnames(dda))
legend("topright", levels(fact$group), fill = pal, bty = "n", xpd = T)


#saving dataset for NA impotation and filtration in NAguideR
write.csv(prot_dda_f, "dda_toNA.csv")
write.csv(prot_dia_f, "dia_toNA.csv")
write.csv(prot_diaml_f, "diaml_toNA.csv")

#We are using NAguideR to find the best way to perform missed values imputation.
#Use generated csv files and "fact_NA.xlsx" as input to NAguidR
library(NAguideR)
NAguideR_app() #In our case, we used CV threshold 0.8, NA ration 0.5 (NA were counted by each group) and used impsqprob as the best method for NA imputation.


#Venn diagram of proteins after filtration (final datasets)
#Further we opening the datasets saved from NAguideR after NA imputation
fin_dia <- data.frame(read.csv("dia_imp.csv", header = TRUE))
fin_diaml <- data.frame(read.csv("diaml_imp.csv", header = TRUE))
fin_dda <- data.frame(read.csv("dda_imp.csv", header = TRUE))

venn.diagram(
  x = list(fin_dia$X, fin_diaml$X, fin_dda$X),
  category.names = c("DIA" , "DIA-ML" , "DDA"),
  filename = '#venn_diagramm_prot_fin.png',
  resolution = 300,
  fill = myCol2,
  print.mode=c("raw","percent"),
  output=F
)


## preparation the transcriptomic data and demonstrating overlaps between proteomics and transcriptomics data
#Opening the data
trans <- data.frame(read.csv("transcriptomics_data.csv", header = TRUE)) # 0 replaced by NA rows with all NA are removed
rownames(trans) <- trans[,1]
trans <- trans[,-1]

#Extracting gene names 
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(trans)
genes <- gsub("\\..*","",genes)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)

#Venn diagrams of proteomics and transcriptomics data overlap
#raw proteins and transcripts
myCol4 <- brewer.pal(4, "Pastel2")
venn.diagram(
  x = list(rownames(prot_dia_f), rownames(prot_diaml_f), rownames(prot_dda_f), G_list$hgnc_symbol),
  category.names = c("DIA" , "DIA-ML" , "DDA", "RNA-seq"),
  filename = '#venn_diagramm_prot_raw_RNAseq.png',
  resolution = 300,
  fill = myCol4,
  print.mode=c("raw","percent"),
  output=F
)

#filtered proteins and transcripts
venn.diagram(
  x = list(fin_dia$X, fin_diaml$X, fin_dda$X, G_list$hgnc_symbol),
  category.names = c("DIA" , "DIA-ML" , "DDA", "RNA-seq"),
  filename = '#venn_diagramm_prot_fin_RNA.png',
  resolution = 300,
  fill = myCol4,
  print.mode=c("raw","percent"),
  output=F
)

#raw peptides and transcripts
venn.diagram(
  x = list(pept_dia$Genes, pept_diaml$Genes, pept_dda$Gene, G_list$hgnc_symbol),
  category.names = c("DIA" , "DIA-ML" , "DDA", "RNA-seq"),
  filename = '#venn_diagramm_pept_genes_RNA-seq.png',
  resolution = 300,
  fill = myCol4,
  print.mode=c("raw","percent"),
  output=F
)


#We using venn function to extract specific unique gene names
library(gplots)
v.table1 <- venn( list(fin_dia$X, fin_diaml$X, fin_dda$X, G_list$hgnc_symbol))
print(v.table1)
