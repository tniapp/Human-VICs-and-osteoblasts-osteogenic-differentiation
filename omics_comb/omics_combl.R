### This part includes analysis of DDA dataset
## Opening the data and performing VSN normalization (normalization was disabled in NAguideR in previous step)
#setwd("your directory")


# OST_VICS_dif
dda_vc_d <- data.frame(read.csv("dda_ost_diffVSvic_diff.csv", header = TRUE)) #upload file after NAguiderR
rownames(dda_vc_d) <- dda_vc_d[,1]

dia_vc_d <- data.frame(read.csv("dia_ost_diffVSvic_diff.csv", header = TRUE)) #upload file after NAguiderR
rownames(dia_vc_d) <- dia_vc_d[,1]

library(dplyr)
common <- intersect(dda_vc_d$X, dia_vc_d$X)  
FCdda <- dda_vc_d[dda_vc_d$X %in% common, c('X', 'logFC')]
FCdda <- FCdda[!duplicated(FCdda$X), ]
FCdda <- arrange(FCdda, FCdda$X)
colnames(FCdda) <- c('Gene_ID', 'logFC_dda')

FCdia <- dia_vc_d[dia_vc_d$X %in% common, c('X', 'logFC')]
FCdia <- FCdia[!duplicated(FCdia$X), ]
FCdia <- arrange(FCdia, FCdia$X)
colnames(FCdia) <- c('Gene_ID', 'logFC_dia')

Dia_dda <- data.frame(FCdda, FCdia)
head(Dia_dda)
library(ggplot2)
#tiff('dda_dia_vics_ost_dif.tiff', units="in", width=16, height=11, res=350, compression = 'lzw')
ggplot(Dia_dda, aes(x=logFC_dda , y=logFC_dia)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold change between DDA and DIA proteomics") + labs(x = "DDA logFC") + labs(y = "DIA logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_dda>1.5 & logFC_dia>1.5,as.character(Gene_ID),'')),hjust=0.5,vjust=-0.5, size = 6) + 
  geom_text(aes(label=ifelse(logFC_dda<(-1.5) & logFC_dia<(-1.5),as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6)
dev.off()

cor(x = Dia_dda$logFC_dda, y = Dia_dda$logFC_dia, method = "pearson")



diaml_vc_d <- data.frame(read.csv("diaml_ost_diffVSvic_diff.csv", header = TRUE)) #upload file after NAguiderR
rownames(diaml_vc_d) <- diaml_vc_d[,1]

common <- intersect(dda_vc_d$X, diaml_vc_d$X)  
FCdda <- dda_vc_d[dda_vc_d$X %in% common, c('X', 'logFC')]
FCdda <- FCdda[!duplicated(FCdda$X), ]
FCdda <- arrange(FCdda, FCdda$X)
colnames(FCdda) <- c('Gene_ID', 'logFC_dda')

FCdiaml <- diaml_vc_d[diaml_vc_d$X %in% common, c('X', 'logFC')]
FCdiaml <- FCdiaml[!duplicated(FCdiaml$X), ]
FCdiaml <- arrange(FCdiaml, FCdiaml$X)
colnames(FCdiaml) <- c('Gene_ID', 'logFC_diaml')

Diaml_dda <- data.frame(FCdda, FCdiaml)
head(Diaml_dda)

#tiff('dda_diaml_vicsVSost_difc.tiff', units="in",width=16, height=11, res=350, compression = 'lzw')
ggplot(Diaml_dda, aes(x=logFC_dda , y=logFC_diaml)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold change between DDA and DIA-ML proteomics") + labs(x = "DDA logFC") + labs(y = "DIA-ML logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_dda>1.5 & logFC_diaml>1.5,as.character(Gene_ID),'')),hjust=0.5,vjust=-0.5, size = 6) + 
  geom_text(aes(label=ifelse(logFC_dda<(-1.5) & logFC_diaml<(-1.5),as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6)
dev.off()

cor(x = Diaml_dda$logFC_dda, y = Diaml_dda$logFC_diaml, method = "pearson")


common <- intersect(dia_vc_d$X, diaml_vc_d$X)  
FCdia <- dia_vc_d[dia_vc_d$X %in% common, c('X', 'logFC')]
FCdia <- FCdia[!duplicated(FCdia$X), ]
FCdia <- arrange(FCdia, FCdia$X)
colnames(FCdia) <- c('Gene_ID', 'logFC_dia')

FCdiaml <- diaml_vc_d[diaml_vc_d$X %in% common, c('X', 'logFC')]
FCdiaml <- FCdiaml[!duplicated(FCdiaml$X), ]
FCdiaml <- arrange(FCdiaml, FCdiaml$X)
colnames(FCdiaml) <- c('Gene_ID', 'logFC_diaml')

Diaml_dia <- data.frame(FCdia, FCdiaml)
head(Diaml_dia)

#tiff('dia_diaml.tiff', units="in", width=16, height=11, res=350, compression = 'lzw')
ggplot(Diaml_dia, aes(x=logFC_dia , y=logFC_diaml)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold change between DIA and DIA-ML proteomics") + labs(x = "DIA logFC") + labs(y = "DIA-ML logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_dia>1.5 & logFC_diaml>1.5,as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6) + 
  geom_text(aes(label=ifelse(logFC_dia<(-1.5) & logFC_diaml<(-1.5),as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6)
dev.off()

cor(x = Diaml_dia$logFC_dia, y = Diaml_dia$logFC_diaml, method = "pearson")



# OST_VICS_contr
dda_vc_c <- data.frame(read.csv("dda_ost_contVSvic_cont.csv", header = TRUE)) #upload file after NAguiderR
rownames(dda_vc_c) <- dda_vc_c[,1]

dia_vc_c <- data.frame(read.csv("dia_ost_contVSvic_cont.csv", header = TRUE)) #upload file after NAguiderR
rownames(dia_vc_c) <- dia_vc_c[,1]

common <- intersect(dda_vc_c$X, dia_vc_c$X)  
FCdda <- dda_vc_c[dda_vc_c$X %in% common, c('X', 'logFC')]
FCdda <- FCdda[!duplicated(FCdda$X), ]
FCdda <- arrange(FCdda, FCdda$X)
colnames(FCdda) <- c('Gene_ID', 'logFC_dda')

FCdia <- dia_vc_c[dia_vc_c$X %in% common, c('X', 'logFC')]
FCdia <- FCdia[!duplicated(FCdia$X), ]
FCdia <- arrange(FCdia, FCdia$X)
colnames(FCdia) <- c('Gene_ID', 'logFC_dia')

Dia_dda <- data.frame(FCdda, FCdia)
head(Dia_dda)

#tiff('dda_dia_vics_ost_contr.tiff', units="in", width=16, height=11, res=350, compression = 'lzw')
ggplot(Dia_dda, aes(x=logFC_dda , y=logFC_dia)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold change between DDA and DIA proteomics ost_vs_vics_contr") + labs(x = "DDA logFC") + labs(y = "DIA logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_dda>1.5 & logFC_dia>1.5,as.character(Gene_ID),'')),hjust=0.5,vjust=-0.5, size = 6) + 
  geom_text(aes(label=ifelse(logFC_dda<(-1.5) & logFC_dia<(-1.5),as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6)
dev.off()

cor(x = Dia_dda$logFC_dda, y = Dia_dda$logFC_dia, method = "pearson")



diaml_vc_c <- data.frame(read.csv("diaml_ost_contVSvic_cont.csv", header = TRUE)) #upload file after NAguiderR
rownames(diaml_vc_c) <- diaml_vc_c[,1]

common <- intersect(dda_vc_c$X, diaml_vc_c$X)  
FCdda <- dda_vc_c[dda_vc_c$X %in% common, c('X', 'logFC')]
FCdda <- FCdda[!duplicated(FCdda$X), ]
FCdda <- arrange(FCdda, FCdda$X)
colnames(FCdda) <- c('Gene_ID', 'logFC_dda')

FCdiaml <- diaml_vc_c[diaml_vc_c$X %in% common, c('X', 'logFC')]
FCdiaml <- FCdiaml[!duplicated(FCdiaml$X), ]
FCdiaml <- arrange(FCdiaml, FCdiaml$X)
colnames(FCdiaml) <- c('Gene_ID', 'logFC_diaml')

Diaml_dda <- data.frame(FCdda, FCdiaml)
head(Diaml_dda)

#tiff('dda_diaml_vicsVSost_contr.tiff', units="in",width=16, height=11, res=350, compression = 'lzw')
ggplot(Diaml_dda, aes(x=logFC_dda , y=logFC_diaml)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold change between DDA and DIA-ML proteomics ostVICscontr") + labs(x = "DDA logFC") + labs(y = "DIA-ML logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_dda>1.5 & logFC_diaml>1.5,as.character(Gene_ID),'')),hjust=0.5,vjust=-0.5, size = 6) + 
  geom_text(aes(label=ifelse(logFC_dda<(-1.5) & logFC_diaml<(-1.5),as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6)
dev.off()

cor(x = Diaml_dda$logFC_dda, y = Diaml_dda$logFC_diaml, method = "pearson")


common <- intersect(dia_vc_c$X, diaml_vc_c$X)  
FCdia <- dia_vc_c[dia_vc_c$X %in% common, c('X', 'logFC')]
FCdia <- FCdia[!duplicated(FCdia$X), ]
FCdia <- arrange(FCdia, FCdia$X)
colnames(FCdia) <- c('Gene_ID', 'logFC_dia')

FCdiaml <- diaml_vc_c[diaml_vc_c$X %in% common, c('X', 'logFC')]
FCdiaml <- FCdiaml[!duplicated(FCdiaml$X), ]
FCdiaml <- arrange(FCdiaml, FCdiaml$X)
colnames(FCdiaml) <- c('Gene_ID', 'logFC_diaml')

Diaml_dia <- data.frame(FCdia, FCdiaml)
head(Diaml_dia)

#tiff('dia_diaml_ostVICs_contr.tiff', units="in", width=16, height=11, res=350, compression = 'lzw')
ggplot(Diaml_dia, aes(x=logFC_dia , y=logFC_diaml)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold change between DIA and DIA-ML proteomics_ostVICs_contr") + labs(x = "DIA logFC") + labs(y = "DIA-ML logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_dia>1.5 & logFC_diaml>1.5,as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6) + 
  geom_text(aes(label=ifelse(logFC_dia<(-1.5) & logFC_diaml<(-1.5),as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6)
dev.off()

cor(x = Diaml_dia$logFC_dia, y = Diaml_dia$logFC_diaml, method = "pearson")



# OST_contr_ost_dif
dda_ost <- data.frame(read.csv("dif_dda_ost_contVSost_diff.csv", header = TRUE)) #upload file after NAguiderR
rownames(dda_ost) <- dda_ost[,1]

dia_ost <- data.frame(read.csv("dif_dia_ost_contVSost_diff.csv", header = TRUE)) #upload file after NAguiderR
rownames(dia_ost) <- dia_ost[,1]

common <- intersect(dda_ost$X, dia_ost$X)  
FCdda <- dda_ost[dda_ost$X %in% common, c('X', 'logFC')]
FCdda <- FCdda[!duplicated(FCdda$X), ]
FCdda <- arrange(FCdda, FCdda$X)
colnames(FCdda) <- c('Gene_ID', 'logFC_dda')

FCdia <- dia_ost[dia_ost$X %in% common, c('X', 'logFC')]
FCdia <- FCdia[!duplicated(FCdia$X), ]
FCdia <- arrange(FCdia, FCdia$X)
colnames(FCdia) <- c('Gene_ID', 'logFC_dia')

Dia_dda <- data.frame(FCdda, FCdia)
head(Dia_dda)

#tiff('dda_dia_ost.tiff', units="in", width=16, height=11, res=350, compression = 'lzw')
ggplot(Dia_dda, aes(x=logFC_dda , y=logFC_dia)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold change between DDA and DIA proteomics ost") + labs(x = "DDA logFC") + labs(y = "DIA logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_dda>1.5 & logFC_dia>1.5,as.character(Gene_ID),'')),hjust=0.5,vjust=-0.5, size = 6) + 
  geom_text(aes(label=ifelse(logFC_dda<(-1.5) & logFC_dia<(-1.5),as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6)
dev.off()

cor(x = Dia_dda$logFC_dda, y = Dia_dda$logFC_dia, method = "pearson")



diaml_ost <- data.frame(read.csv("dif_diaml_ost_contVSost_diff.csv", header = TRUE)) #upload file after NAguiderR
rownames(diaml_ost) <- diaml_ost[,1]

common <- intersect(dda_ost$X, diaml_ost$X)  
FCdda <- dda_ost[dda_ost$X %in% common, c('X', 'logFC')]
FCdda <- FCdda[!duplicated(FCdda$X), ]
FCdda <- arrange(FCdda, FCdda$X)
colnames(FCdda) <- c('Gene_ID', 'logFC_dda')

FCdiaml <- diaml_ost[diaml_ost$X %in% common, c('X', 'logFC')]
FCdiaml <- FCdiaml[!duplicated(FCdiaml$X), ]
FCdiaml <- arrange(FCdiaml, FCdiaml$X)
colnames(FCdiaml) <- c('Gene_ID', 'logFC_diaml')

Diaml_dda <- data.frame(FCdda, FCdiaml)
head(Diaml_dda)

#tiff('dda_diaml_ost.tiff', units="in",width=16, height=11, res=350, compression = 'lzw')
ggplot(Diaml_dda, aes(x=logFC_dda , y=logFC_diaml)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold change between DDA and DIA-ML proteomics osts") + labs(x = "DDA logFC") + labs(y = "DIA-ML logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_dda>1.5 & logFC_diaml>1.5,as.character(Gene_ID),'')),hjust=0.5,vjust=-0.5, size = 6) + 
  geom_text(aes(label=ifelse(logFC_dda<(-1.5) & logFC_diaml<(-1.5),as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6)
dev.off()

cor(x = Diaml_dda$logFC_dda, y = Diaml_dda$logFC_diaml, method = "pearson")


common <- intersect(dia_ost$X, diaml_ost$X)  
FCdia <- dia_ost[dia_ost$X %in% common, c('X', 'logFC')]
FCdia <- FCdia[!duplicated(FCdia$X), ]
FCdia <- arrange(FCdia, FCdia$X)
colnames(FCdia) <- c('Gene_ID', 'logFC_dia')

FCdiaml <- diaml_ost[diaml_ost$X %in% common, c('X', 'logFC')]
FCdiaml <- FCdiaml[!duplicated(FCdiaml$X), ]
FCdiaml <- arrange(FCdiaml, FCdiaml$X)
colnames(FCdiaml) <- c('Gene_ID', 'logFC_diaml')

Diaml_dia <- data.frame(FCdia, FCdiaml)
head(Diaml_dia)

#tiff('dia_diaml_ost.tiff', units="in", width=16, height=11, res=350, compression = 'lzw')
ggplot(Diaml_dia, aes(x=logFC_dia , y=logFC_diaml)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold change between DIA and DIA-ML proteomics_ost") + labs(x = "DIA logFC") + labs(y = "DIA-ML logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_dia>1.5 & logFC_diaml>1.5,as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6) + 
  geom_text(aes(label=ifelse(logFC_dia<(-1.5) & logFC_diaml<(-1.5),as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6)
dev.off()

cor(x = Diaml_dia$logFC_dia, y = Diaml_dia$logFC_diaml, method = "pearson")


# VICS_contr_VICs_dif
dda_vic <- data.frame(read.csv("dda_vic_contVSvic_diff.csv", header = TRUE)) #upload file after NAguiderR
rownames(dda_vic) <- dda_vic[,1]

dia_vic <- data.frame(read.csv("dia_vic_contVSvic_diff.csv", header = TRUE)) #upload file after NAguiderR
rownames(dia_vic) <- dia_vic[,1]

common <- intersect(dda_vic$X, dia_vic$X)  
FCdda <- dda_vic[dda_vic$X %in% common, c('X', 'logFC')]
FCdda <- FCdda[!duplicated(FCdda$X), ]
FCdda <- arrange(FCdda, FCdda$X)
colnames(FCdda) <- c('Gene_ID', 'logFC_dda')

FCdia <- dia_vic[dia_vic$X %in% common, c('X', 'logFC')]
FCdia <- FCdia[!duplicated(FCdia$X), ]
FCdia <- arrange(FCdia, FCdia$X)
colnames(FCdia) <- c('Gene_ID', 'logFC_dia')

Dia_dda <- data.frame(FCdda, FCdia)
head(Dia_dda)

#tiff('dda_dia_vic.tiff', units="in", width=16, height=11, res=350, compression = 'lzw')
ggplot(Dia_dda, aes(x=logFC_dda , y=logFC_dia)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold change between DDA and DIA proteomics vic") + labs(x = "DDA logFC") + labs(y = "DIA logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_dda>1.5 & logFC_dia>1.5,as.character(Gene_ID),'')),hjust=0.5,vjust=-0.5, size = 6) + 
  geom_text(aes(label=ifelse(logFC_dda<(-1.5) & logFC_dia<(-1.5),as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6)
dev.off()

cor(x = Dia_dda$logFC_dda, y = Dia_dda$logFC_dia, method = "pearson")



diaml_vic <- data.frame(read.csv("diaml_vic_contVSvic_diff.csv", header = TRUE)) #upload file after NAguiderR
rownames(diaml_vic) <- diaml_vic[,1]

common <- intersect(dda_vic$X, diaml_vic$X)  
FCdda <- dda_vic[dda_vic$X %in% common, c('X', 'logFC')]
FCdda <- FCdda[!duplicated(FCdda$X), ]
FCdda <- arrange(FCdda, FCdda$X)
colnames(FCdda) <- c('Gene_ID', 'logFC_dda')

FCdiaml <- diaml_vic[diaml_vic$X %in% common, c('X', 'logFC')]
FCdiaml <- FCdiaml[!duplicated(FCdiaml$X), ]
FCdiaml <- arrange(FCdiaml, FCdiaml$X)
colnames(FCdiaml) <- c('Gene_ID', 'logFC_diaml')

Diaml_dda <- data.frame(FCdda, FCdiaml)
head(Diaml_dda)

#tiff('dda_diaml_vic.tiff', units="in",width=16, height=11, res=350, compression = 'lzw')
ggplot(Diaml_dda, aes(x=logFC_dda , y=logFC_diaml)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold change between DDA and DIA-ML proteomics vic") + labs(x = "DDA logFC") + labs(y = "DIA-ML logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_dda>1.5 & logFC_diaml>1.5,as.character(Gene_ID),'')),hjust=0.5,vjust=-0.5, size = 6) + 
  geom_text(aes(label=ifelse(logFC_dda<(-1.5) & logFC_diaml<(-1.5),as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6)
dev.off()

cor(x = Diaml_dda$logFC_dda, y = Diaml_dda$logFC_diaml, method = "pearson")


common <- intersect(dia_vic$X, diaml_vic$X)  
FCdia <- dia_vic[dia_vic$X %in% common, c('X', 'logFC')]
FCdia <- FCdia[!duplicated(FCdia$X), ]
FCdia <- arrange(FCdia, FCdia$X)
colnames(FCdia) <- c('Gene_ID', 'logFC_dia')

FCdiaml <- diaml_vic[diaml_vic$X %in% common, c('X', 'logFC')]
FCdiaml <- FCdiaml[!duplicated(FCdiaml$X), ]
FCdiaml <- arrange(FCdiaml, FCdiaml$X)
colnames(FCdiaml) <- c('Gene_ID', 'logFC_diaml')

Diaml_dia <- data.frame(FCdia, FCdiaml)
head(Diaml_dia)

#tiff('dia_diaml_vic.tiff', units="in", width=16, height=11, res=350, compression = 'lzw')
ggplot(Diaml_dia, aes(x=logFC_dia , y=logFC_diaml)) +
  geom_point(shape= 16) + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  labs(title = "Comparison of gene fold change between DIA and DIA-ML proteomics_vic") + labs(x = "DIA logFC") + labs(y = "DIA-ML logFC") +
  theme(legend.position = "top") +
  geom_text(aes(label=ifelse(logFC_dia>1.5 & logFC_diaml>1.5,as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6) + 
  geom_text(aes(label=ifelse(logFC_dia<(-1.5) & logFC_diaml<(-1.5),as.character(Gene_ID),'')),hjust=0,vjust=-0.5, size = 6)
dev.off()

cor(x = Diaml_dia$logFC_dia, y = Diaml_dia$logFC_diaml, method = "pearson")







#Venn diagramm
up.genes_dda <- dda_vc_d[dda_vc_d$logFC > 1 & dda_vc_d$adj.P.Val < 0.05, 1] 
up.genes_dia <- dia_vc_d[dia_vc_d$logFC > 1 & dia_vc_d$adj.P.Val < 0.05, 1] 
up.genes_dia_ml <- diaml_vc_d[diaml_vc_d$logFC > 1 & diaml_vc_d$adj.P.Val < 0.05, 1] 

dn.genes_dda <- dda_vc_d[dda_vc_d$logFC < 1 & dda_vc_d$adj.P.Val < 0.05, 1] 
dn.genes_dia <- dia_vc_d[dia_vc_d$logFC < 1 & dia_vc_d$adj.P.Val < 0.05, 1] 
dn.genes_dia_ml <- diaml_vc_d[diaml_vc_d$logFC < 1 & diaml_vc_d$adj.P.Val < 0.05, 1] 

library(VennDiagram)
library(RColorBrewer)
myCol1 <- brewer.pal(3, "Pastel2")

#This code will draw diagram to working directory
venn.diagram(
  x = list(up.genes_dda, up.genes_dia, up.genes_dia_ml), 
  category.names = c("DDA" , "DIA" , "DIA-ML"),
  filename = '#venn_diagramm_up.png',
  resolution = 600,
  fill = myCol1,
  output=F
)

venn.diagram(
  x = list(dn.genes_dda, dn.genes_dia, dn.genes_dia_ml), 
  category.names = c("DDA" , "DIA" , "DIA-ML"),
  filename = '#venn_diagramm_dn.png',
  resolution = 600,
  fill = myCol1,
  output=F
)

