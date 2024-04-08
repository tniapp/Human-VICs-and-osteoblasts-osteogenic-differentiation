### This part includes analysis of DDA dataset
## Opening the data and performing VSN normalization (normalization was disabled in NAguideR in previous step)
#setwd("your directory")
dia <- data.frame(read.csv("dia_imp.csv", header = TRUE)) #upload file after NAguiderR
rownames(dia) <- dia[,1]
dia <- dia[,-1]

#create sample information file:

samples <- colnames(dia)
samples

group <- c('ost_cont','ost_diff','ost_cont','ost_diff','ost_cont', 'ost_diff',
          'ost_cont','ost_diff','ost_cont','ost_diff','ost_cont', 'ost_diff',
          'vic_cont','vic_diff','vic_cont','vic_diff','vic_cont',
          'vic_diff','vic_cont','vic_diff','vic_cont','vic_diff')

differentiation <- c('cont','diff','cont','diff','cont','diff',
                    'cont','diff','cont','diff','cont','diff',
                    'cont','diff','cont','diff','cont',
                    'diff','cont','diff','cont','diff') 

cell_type <- c('ost','ost','ost','ost','ost','ost','ost','ost','ost','ost','ost','ost',
         'vic','vic','vic','vic','vic','vic','vic','vic','vic','vic')

donor <- c('8','8','9','9','11','11','12','12','20','20','22','22',
           '213','213','431','431','432','432','482','482','494','494')


fact <- data.frame(samples,differentiation,cell_type,donor,group)

fact$differentiation <- as.factor(fact$differentiation)
fact$differentiatio
fact$cell_type <- as.factor(fact$cell_type)
fact$cell_type
fact$donor <- as.factor(fact$donor)
fact$donor
fact$group <- as.factor(fact$group)
fact$group

rownames(fact) <- fact[,1]
fact

##Data normalization
library(vsn)
library(RColorBrewer)
pal <- brewer.pal(n = 4, name = "Set1")
cols <- pal[fact$group]
#Log-transformed data without normalization
boxplot(dia, outline = TRUE, col = cols, main = "DIA data intensity distribution after NAguideR", names=colnames(dia))
legend("topright", levels(fact$group), fill = pal, bty = "n", xpd = T)

maplot <- function(X1, X2, pch = 21, main = "MA-plot", xlab = "Average log-expression", ylab = "Expression log-ratio", lpars = list(col = "blue", lwd = 2), ...){
  X <- (rowMeans(X2) + rowMeans(X1)) / 2
  Y <- rowMeans(X2) - rowMeans(X1)
  scatter.smooth(x = X, y = Y,
                 main = main, pch = pch,
                 xlab = xlab, ylab = ylab,
                 lpars = lpars, ...)
  abline(h = c(-1, 0, 1), lty = c(2, 1, 2))
}

maplot(dia[, fact$differentiation == "cont"], dia[, fact$differentiation == "diff"], main = "Log-expression data")
meanSdPlot(as.matrix(dia))

#VSN normalization
library(limma)
dia_vsn <- normalizeVSN(dia)

boxplot(dia_vsn, outline = TRUE, col = cols, main = "DIA data intensity distribution after VSN", names=colnames(dia))
legend("topright", levels(fact$group), fill = pal, bty = "n", xpd = T)
maplot(dia_vsn[, fact$differentiation == "cont"], dia_vsn[, fact$differentiation == "diff"], main = "Log-expression data")
meanSdPlot(as.matrix(dia_vsn))

#quantile normalization
dia_Q <- normalizeQuantiles(dia)

boxplot(dia_Q, outline = TRUE, col = cols, main = "DIA data intensity distribution after normalizeQuantiles", names=colnames(dia))
legend("topright", levels(fact$group), fill = pal, bty = "n", xpd = T)
maplot(dia_Q[, fact$differentiation == "cont"], dia_Q[, fact$differentiation == "diff"], main = "Log-expression data")
meanSdPlot(as.matrix(dia_Q))

#Quantile noemalization seems more accurate, so we will use it further
dia <- dia_Q

## multivariate analysis
#nMDS
library(vegan)
set.seed(125)
nmds <- metaMDS(t(dia), distance = "euclidean")
nmds

data.scores <- as.data.frame(scores(nmds)$sites)

data.scores$group <- fact$group
head(data.scores)

library(ggplot2)

#tiff('dia_nMDS.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = group))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "group", y = "NMDS2") + scale_colour_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"))
#dev.off()

#PCA
library(mixOmics)

dia_pca <- pca(t(dia), ncomp = 8)
plot(dia_pca)

dia_pca <- pca(t(dia), ncomp = 3)
#tiff('dia_pca.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
plotIndiv(dia_pca, comp = c(1,2), legend = TRUE,
          group = fact$group,
          ellipse = TRUE,
          title = 'DIA PCA comp 1 - 2')
#dev.off()


#tiff('dda_pca_donor.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
plotIndiv(dia_pca, comp = c(1,2), legend = TRUE,
          group = fact$donor,
          ellipse = TRUE,
          title = 'DDA PCA comp 1 - 2')
#dev.off()


#PLS-DA
ordination.optimum.splsda <- splsda(t(dia), fact$group, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#tiff('dia_splsda.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 1.5, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, title = "PLS-DA ordination", size.title = 1.5, legend=TRUE)
#dev.off()
layout(1,1)


##Differential expression analysis (see guide https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf)
design = model.matrix(~0+fact$group)
colnames(design) = c('ost_cont', 'ost_diff','vic_cont','vic_diff')

#To make all pair-wise comparisons between all four groups one could proceed

fit_DIA <- lmFit(dia, design)
contrasts_groups = c('ost_cont-ost_diff','ost_cont-vic_cont','vic_cont-vic_diff','ost_diff-vic_diff')
contrast.matrix <- makeContrasts(contrasts = contrasts_groups,levels=design)

fit2_DIA <- contrasts.fit(fit_DIA, contrast.matrix)
fit2_DIA <- eBayes(fit2_DIA)

topTable(fit2_DIA, coef=1, adjust="BH")
results_DIA <- decideTests(fit2_DIA)

#tiff('venn_all_comp_DIA.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
vennDiagram(results_DIA)
#dev.off()

#results save
library(EnhancedVolcano)

dia_ost_contVSost_diff <- topTable(fit2_DIA, coef=1, adjust="BH", number = as.numeric(length(dia$Ost_contr_11)))
#write.csv(dia_ost_contVSost_diff, "dif_dia_ost_contVSost_diff.csv")

#tiff('dia_ost_contVSost_diff.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
EnhancedVolcano(dia_ost_contVSost_diff,
                lab = rownames(dia_ost_contVSost_diff),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
               xlim = c(-2.5, 2), # needs to be adjusted annually!
               ylim = c(0, 5.2), # needs to be adjusted annually!
                FCcutoff = 1,
                title ="DIA volcano plot for 'ost_cont-ost_diff'",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()


dia_ost_contVSvic_cont <- topTable(fit2_DIA, coef=2, adjust="BH", number = as.numeric(length(dia$Ost_contr_11)))
#write.csv(dia_ost_contVSvic_cont, "dia_ost_contVSvic_cont.csv")

#tiff('dia_ost_contVSvic_cont.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
EnhancedVolcano(dia_ost_contVSvic_cont,
                lab = rownames(dia_ost_contVSvic_cont),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-5, 3.5), # needs to be adjusted annually!
                ylim = c(0, 10), # needs to be adjusted annually!
                FCcutoff = 1,
                title ="DIA volcano plot for 'dda_ost_cont-vic_cont'",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()


dia_vic_contVSvic_diff <- topTable(fit2_DIA, coef=3, adjust="BH", number = as.numeric(length(dia$Ost_contr_11)))
#write.csv(dia_vic_contVSvic_diff, "dia_vic_contVSvic_diff.csv")

#tiff('dia_vic_contVSvic_diff.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
EnhancedVolcano(dia_vic_contVSvic_diff,
                lab = rownames(dia_vic_contVSvic_diff),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-2.5, 2.2), # needs to be adjusted annually!
                ylim = c(0, 3.5), # needs to be adjusted annually!
                FCcutoff = 1,
                title ="DIA volcano plot for 'vic_cont-vic_diff",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()


dia_ost_diffVSvic_diff <- topTable(fit2_DIA, coef=4, adjust="BH", number = as.numeric(length(dia$Ost_contr_11)))
#write.csv(dia_ost_diffVSvic_diff, "dia_ost_diffVSvic_diff.csv")

#tiff('dia_ost_diffVSvic_diff.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
EnhancedVolcano(dia_ost_diffVSvic_diff,
                lab = rownames(dia_ost_diffVSvic_diff),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
              xlim = c(-5, 4), # needs to be adjusted annually!
               ylim = c(0, 11), # needs to be adjusted annually!
                FCcutoff = 1,
                title ="DIA volcano plot for 'dda_ost_diff-vic_diff'",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()



diam <- as.matrix(dia)

###!!! not yet done
#Top-proteins verification
#tiff('ost_vic_sim_dia.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
par(mfrow=c(2,2))
boxplot(diam[c("SAMHD1"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1) #not found
boxplot(diam[c("MAOA"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1) #not found
boxplot(diam[c("METTL7A"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1) #not found
boxplot(diam[c("ALCAM"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
#dev.off()

#tiff('ost_vic_dif_dia_1.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
par(mfrow=c(2,2))
boxplot(diam[c("COL1A1"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1) #not found
boxplot(diam[c("MYH10"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1) #not found
boxplot(diam[c("HTRA1"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(diam[c("ITGA2"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
#dev.off()

#tiff('ost_vic_dif_dia_2.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
par(mfrow=c(2,2))
boxplot(diam[c("TGM2"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1) #not found
boxplot(diam[c("TGFBI"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1) #not found
boxplot(diam[c("STAT2"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1) 
boxplot(diam[c("METTL7A"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1) #not found
#dev.off()

layout(1,1)

