### This part includes analysis of DDA dataset
## Opening the data and performing VSN normalization (normalization was disabled in NAguideR in previous step)
#setwd("your directory")
dda <- data.frame(read.csv("dda_imp.csv", header = TRUE)) #upload file after NAguiderR
rownames(dda) <- dda[,1]
dda <- dda[,-1]

#create sample information file:

samples <- colnames(dda)
samples

group <- c('ost_cont','ost_cont','ost_cont','ost_cont','ost_cont', 'ost_cont',
          'ost_diff','ost_diff','ost_diff','ost_diff','ost_diff', 'ost_diff',
          'vic_cont','vic_cont','vic_cont','vic_cont','vic_cont',
          'vic_diff','vic_diff','vic_diff','vic_diff','vic_diff')

differentiation <- c('cont','cont','cont','cont','cont','cont',
                    'diff','diff','diff','diff','diff','diff',
                    'cont','cont','cont','cont','cont',
                    'diff','diff','diff','diff','diff') 

cell_type <- c('ost','ost','ost','ost','ost','ost','ost','ost','ost','ost','ost','ost',
         'vic','vic','vic','vic','vic','vic','vic','vic','vic','vic')

donor <- c('11','12','20','22','8','9','11','12','20','22','8','9',
           '213','431','432','482','494','213','431','432','482','494')


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
boxplot(dda, outline = TRUE, col = cols, main = "DDA data intensity distribution after NAguideR", names=colnames(dda))
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

maplot(dda[, fact$differentiation == "cont"], dda[, fact$differentiation == "diff"], main = "Log-expression data")
meanSdPlot(as.matrix(dda))

#VSN normalization
library(limma)
dda_vsn <- normalizeVSN(dda)

boxplot(dda_vsn, outline = TRUE, col = cols, main = "DDA data intensity distribution after VSN", names=colnames(dda))
legend("topright", levels(fact$group), fill = pal, bty = "n", xpd = T)
maplot(dda_vsn[, fact$differentiation == "cont"], dda_vsn[, fact$differentiation == "diff"], main = "Log-expression data")
meanSdPlot(as.matrix(dda_vsn))

#quantile normalization
dda_Q <- normalizeQuantiles(dda)

boxplot(dda_Q, outline = TRUE, col = cols, main = "DDA data intensity distribution after normalizeQuantiles", names=colnames(dda))
legend("topright", levels(fact$group), fill = pal, bty = "n", xpd = T)
maplot(dda_Q[, fact$differentiation == "cont"], dda_Q[, fact$differentiation == "diff"], main = "Log-expression data")
meanSdPlot(as.matrix(dda_Q))

#Quantile noemalization seems more accurate, so we will use it further
dda <- dda_Q

## multivariate analysis
#nMDS
library(vegan)
set.seed(125)
nmds <- metaMDS(t(dda), distance = "euclidean")
nmds

data.scores <- as.data.frame(scores(nmds)$sites)

data.scores$group <- fact$group
head(data.scores)

library(ggplot2)
#tiff('dda_nMDS.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
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

dda_pca <- pca(t(dda), ncomp = 8)
plot(dda_pca)

dda_pca <- pca(t(dda), ncomp = 3)
#tiff('dda_pca.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
plotIndiv(dda_pca, comp = c(1,2), legend = TRUE,
          group = fact$group,
          ellipse = TRUE,
          title = 'DDA PCA comp 1 - 2')
#dev.off()


#tiff('dda_pca_donor.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
plotIndiv(dda_pca, comp = c(1,2), legend = TRUE,
          group = fact$donor,
          ellipse = TRUE,
          title = 'DDA PCA comp 1 - 2')
#dev.off()


#PLS-DA
#PLS-DA 
ordination.optimum.splsda <- splsda(t(dda), fact$group, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#tiff('dda_splsda.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
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

fit_DDA <- lmFit(dda, design)
contrasts_groups = c('ost_cont-ost_diff','ost_cont-vic_cont','vic_cont-vic_diff','ost_diff-vic_diff')
contrast.matrix <- makeContrasts(contrasts = contrasts_groups,levels=design)

fit2_DDA <- contrasts.fit(fit_DDA, contrast.matrix)
fit2_DDA <- eBayes(fit2_DDA)

topTable(fit2_DDA, coef=1, adjust="BH")
results_DDA <- decideTests(fit2_DDA)

#tiff('venn_all_comp.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
vennDiagram(results_DDA)
#dev.off()

#results save
library(EnhancedVolcano)

dda_ost_contVSost_diff <- topTable(fit2_DDA, coef=1, adjust="BH", number = as.numeric(length(dda$Ost_contr_11)))
#write.csv(dda_ost_contVSost_diff, "dif_dda_ost_contVSost_diff.csv")

#tiff('dda_ost_contVSost_diff.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
EnhancedVolcano(dda_ost_contVSost_diff,
                lab = rownames(dda_ost_contVSost_diff),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
               xlim = c(-3, 2), # needs to be adjusted annually!
               ylim = c(0, 8), # needs to be adjusted annually!
                FCcutoff = 1,
                title ="DDA volcano plot for 'ost_cont-ost_diff'",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()


dda_ost_contVSvic_cont <- topTable(fit2_DDA, coef=2, adjust="BH", number = as.numeric(length(dda$Ost_contr_11)))
#write.csv(dda_ost_contVSvic_cont, "dda_ost_contVSvic_cont.csv")

#tiff('dda_ost_contVSvic_cont.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
EnhancedVolcano(dda_ost_contVSvic_cont,
                lab = rownames(dda_ost_contVSvic_cont),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-4, 2.8), # needs to be adjusted annually!
                ylim = c(0, 8), # needs to be adjusted annually!
                FCcutoff = 1,
                title ="DDA volcano plot for 'dda_ost_cont-vic_cont'",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()


dda_vic_contVSvic_diff <- topTable(fit2_DDA, coef=3, adjust="BH", number = as.numeric(length(dda$Ost_contr_11)))
#write.csv(dda_vic_contVSvic_diff, "dda_vic_contVSvic_diff.csv")

#tiff('dda_vic_contVSvic_diff.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
EnhancedVolcano(dda_vic_contVSvic_diff,
                lab = rownames(dda_vic_contVSvic_diff),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-3, 2.5), # needs to be adjusted annually!
                ylim = c(0, 7.5), # needs to be adjusted annually!
                FCcutoff = 1,
                title ="DDA volcano plot for 'vic_cont-vic_diff",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
# dev.off()


dda_ost_diffVSvic_diff <- topTable(fit2_DDA, coef=4, adjust="BH", number = as.numeric(length(dda$Ost_contr_11)))
#write.csv(dda_ost_diffVSvic_diff, "dda_ost_diffVSvic_diff.csv")

#tiff('dda_ost_diffVSvic_diff.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
EnhancedVolcano(dda_ost_diffVSvic_diff,
                lab = rownames(dda_ost_diffVSvic_diff),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
               xlim = c(-4.5, 2.5), # needs to be adjusted annually!
                ylim = c(0, 10), # needs to be adjusted annually!
                FCcutoff = 1,
                title ="DDA volcano plot for 'dda_ost_diff-vic_diff'",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()



ddam <- as.matrix(dda)

#Top-proteins verification
#tiff('ost_vic_sim_dda.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
par(mfrow=c(2,2))
boxplot(ddam[c("SAMHD1"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(ddam[c("MAOA"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(ddam[c("METTL7A"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(ddam[c("ALCAM"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
#dev.off()

#tiff('ost_vic_dif_dda_1.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
par(mfrow=c(2,2))
boxplot(ddam[c("COL1A1"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(ddam[c("MYH10"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(ddam[c("HTRA1"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(ddam[c("ITGA2"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
#dev.off()

#tiff('ost_vic_dif_dda_2.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
par(mfrow=c(2,2))
boxplot(ddam[c("TGM2"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(ddam[c("TGFBI"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(ddam[c("STAT2"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(ddam[c("METTL7A"),] ~ group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
#dev.off()

layout(1,1)

