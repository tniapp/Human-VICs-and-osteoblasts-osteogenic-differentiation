
dda <- data.frame(read.csv("dda_ost_contVSvic_cont.csv", header = TRUE)) 
dia <- data.frame(read.csv("dia_ost_contVSvic_cont.csv", header = TRUE)) 
diaml <- data.frame(read.csv("diaml_ost_contVSvic_cont.csv", header = TRUE)) 

dda_up <- dda[dda$logFC > 1 & dda$adj.P.Val < 0.05, 1] 
dda_dn <- dda[dda$logFC < -1 & dda$adj.P.Val < 0.05, 1] 

dia_up <- dia[dia$logFC > 1 & dia$adj.P.Val < 0.05, 1] 
dia_dn <- dia[dia$logFC < -1 & dia$adj.P.Val < 0.05, 1] 

diaml_up <- diaml[diaml$logFC > 1 & diaml$adj.P.Val < 0.05, 1] 
diaml_dn <- diaml[diaml$logFC < -1 & diaml$adj.P.Val < 0.05, 1] 



#Enrichment plot for upregulated
ProtW_up <- c(dia_up, diaml_up, dda_up)
ProtW_up <- unique(ProtW_up)
ProtRel_up <- VennDiagram::get.venn.partitions(list(dia=dia_up, diaml=diaml_up, dda=dda_up))
ProtRel_up <- ProtRel_up[1,5]
ProtRel_up <- as.data.frame(ProtRel_up)
ProtRel_up <- ProtRel_up$X1

library("org.Hs.eg.db")
library(rWikiPathways)
library(clusterProfiler)

ProtW_up <- clusterProfiler::bitr(ProtW_up,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
ProtRel_up <- clusterProfiler::bitr(ProtRel_up,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


egobp_ProtW_up <- clusterProfiler::enrichGO(
  gene     = ProtW_up[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05, 
  readable = TRUE)

egobp_ProtRel_up <- clusterProfiler::enrichGO(
  gene     = ProtRel_up[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE)


dotplot(egobp_ProtW_up, showCategory = 20)
dotplot(egobp_ProtRel_up, showCategory = 20)

egobp_ProtW_up <- data.frame(egobp_ProtW_up)
egobp_ProtRel_up <- data.frame(egobp_ProtRel_up)

egobp_ProtW_up <- egobp_ProtW_up[,c(2,6,9)]
egobp_ProtW_up$Dataset <- "Prot. w."
egobp_ProtRel_up <- egobp_ProtRel_up[,c(2,6,9)]
egobp_ProtRel_up$Dataset <- "Prot. rel."
egobp <- rbind(egobp_ProtW_up, egobp_ProtRel_up)
head(egobp)

library(viridis)
library(hrbrthemes)
library(ggplot2)

#tiff('ost_vic_contr_up_bp.tiff', units="in", width=7, height=7, res=300, compression = 'lzw')
ggplot(egobp, aes(x= Dataset, y=Description, size=Count, color=p.adjust)) + 
  geom_point(alpha = 0.8) +
  scale_color_viridis() +
  theme_bw()
dev.off()

#Enrichment plot for downregulated
ProtW_dn <- c(dia_dn, diaml_dn, dda_dn)
ProtW_dn <- unique(ProtW_dn)
ProtRel_dn <- VennDiagram::get.venn.partitions(list(dia=dia_dn, diaml=diaml_dn, dda=dda_dn))
ProtRel_dn <- ProtRel_dn[1,5]
ProtRel_dn <- as.data.frame(ProtRel_dn)
ProtRel_dn <- ProtRel_dn$X1

ProtW_dn <- clusterProfiler::bitr(ProtW_dn,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
ProtRel_dn <- clusterProfiler::bitr(ProtRel_dn,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


egobp_ProtW_dn <- clusterProfiler::enrichGO(
  gene     = ProtW_dn[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05, 
  readable = TRUE)

egobp_ProtRel_dn <- clusterProfiler::enrichGO(
  gene     = ProtRel_dn[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE)


egobp_ProtW_dn <- data.frame(egobp_ProtW_dn)
egobp_ProtRel_dn <- data.frame(egobp_ProtRel_dn)
head(egobp_ProtW_dn)


egobp_ProtW_dn <- egobp_ProtW_dn[,c(2,6,9)]
egobp_ProtW_dn$Dataset <- "Prot. w."
egobp_ProtRel_dn <- egobp_ProtRel_dn[,c(2,6,9)]
egobp_ProtRel_dn$Dataset <- "Prot. rel."
egobp_dn <- rbind(egobp_ProtW_dn, egobp_ProtRel_dn)
head(egobp_dn)


#tiff('ost_vic_contr_dn_bp.tiff', units="in", width=7, height=7, res=300, compression = 'lzw')
ggplot(egobp_dn, aes(x= Dataset, y=Description, size=Count, color=p.adjust)) + 
  geom_point(alpha = 0.8) +
  scale_color_viridis() +
  theme_bw()
dev.off()



#########DIFFERENTIATED cells

dda_dif <- data.frame(read.csv("dda_ost_diffVSvic_diff.csv", header = TRUE)) 
dia_dif <- data.frame(read.csv("dia_ost_diffVSvic_diff.csv", header = TRUE)) 
diaml_dif <- data.frame(read.csv("diaml_ost_diffVSvic_diff.csv", header = TRUE)) 
dda_dif_up <- dda_dif[dda_dif$logFC > 1 & dda_dif$adj.P.Val < 0.05, 1] 
dda_dif_dn <- dda_dif[dda_dif$logFC < -1 & dda_dif$adj.P.Val < 0.05, 1] 

dia_dif_up <- dia_dif[dia_dif$logFC > 1 & dia_dif$adj.P.Val < 0.05, 1] 
dia_dif_dn <- dia_dif[dia_dif$logFC < -1 & dia_dif$adj.P.Val < 0.05, 1] 

diaml_dif_up <- diaml_dif[diaml_dif$logFC > 1 & diaml_dif$adj.P.Val < 0.05, 1] 
diaml_dif_dn <- diaml_dif[diaml_dif$logFC < -1 & diaml_dif$adj.P.Val < 0.05, 1] 



#Enrichment plot for upregulated
ProtW_dif_up <- c(dia_dif_up, diaml_dif_up, dda_dif_up)
ProtW_dif_up <- unique(ProtW_dif_up)
ProtRel_dif_up <- VennDiagram::get.venn.partitions(list(dia=dia_dif_up, diaml=diaml_dif_up, dda=dda_dif_up))
ProtRel_dif_up <- ProtRel_dif_up[1,5]
ProtRel_dif_up <- as.data.frame(ProtRel_dif_up)
ProtRel_dif_up <- ProtRel_dif_up$X1

ProtW_dif_up <- clusterProfiler::bitr(ProtW_dif_up,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
ProtRel_dif_up <- clusterProfiler::bitr(ProtRel_dif_up,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


egobp_ProtW_dif_up <- clusterProfiler::enrichGO(
  gene     = ProtW_dif_up[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05, 
  readable = TRUE)

egobp_ProtRel_dif_up <- clusterProfiler::enrichGO(
  gene     = ProtRel_dif_up[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE)


dotplot(egobp_ProtW_dif_up, showCategory = 20)
dotplot(egobp_ProtRel_dif_up, showCategory = 20)

egobp_ProtW_dif_up <- data.frame(egobp_ProtW_dif_up)
egobp_ProtRel_dif_up <- data.frame(egobp_ProtRel_dif_up)

egobp_ProtW_dif_up <- egobp_ProtW_dif_up[,c(2,6,9)]
egobp_ProtW_dif_up$Dataset <- "Prot. w.dif"
egobp_ProtRel_dif_up <- egobp_ProtRel_dif_up[,c(2,6,9)]
egobp_ProtRel_dif_up$Dataset <- "Prot. rel.dif"
egobp_dif <- rbind(egobp_ProtW_dif_up, egobp_ProtRel_dif_up)
head(egobp_dif)


#tiff('ost_vic_dif_up_bp.tiff', units="in", width=10, height=12, res=300, compression = 'lzw')
ggplot(egobp_dif, aes(x= Dataset, y=Description, size=Count, color=p.adjust)) + 
  geom_point(alpha = 0.8) +
  scale_color_viridis() +
  theme_bw()
dev.off()

#Enrichment plot for downregulated
ProtW_dif_dn <- c(dia_dif_dn, diaml_dif_dn, dda_dif_dn)
ProtW_dif_dn <- unique(ProtW_dif_dn)
ProtRel_dif_dn <- VennDiagram::get.venn.partitions(list(dia=dia_dif_dn, diaml=diaml_dif_dn, dda=dda_dif_dn))
ProtRel_dif_dn <- ProtRel_dif_dn[1,5]
ProtRel_dif_dn <- as.data.frame(ProtRel_dif_dn)
ProtRel_dif_dn <- ProtRel_dif_dn$X1

ProtW_dif_dn <- clusterProfiler::bitr(ProtW_dif_dn,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
ProtRel_dif_dn <- clusterProfiler::bitr(ProtRel_dif_dn,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


egobp_ProtW_dif_dn <- clusterProfiler::enrichGO(
  gene     = ProtW_dif_dn[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05, 
  readable = TRUE)

egobp_ProtRel_dif_dn <- clusterProfiler::enrichGO(
  gene     = ProtRel_dif_dn[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE)


egobp_ProtW_dif_dn <- data.frame(egobp_ProtW_dif_dn)
egobp_ProtRel_dif_dn <- data.frame(egobp_ProtRel_dif_dn)

egobp_ProtRel_dif_dn <- egobp_ProtRel_dif_dn[,c(1,2,3)]
egobp_ProtRel_dif_dn$Dataset <- "ProtRel_dif_dn"
egobp_ProtW_dif_dn <- egobp_ProtW_dif_dn[,c(2,6,9)]
egobp_ProtW_dif_dn$Dataset <- "ProtW_dif_dn"

egobp_dif_dn <- rbind(egobp_ProtRel_dif_dn, egobp_ProtW_dif_dn)
head(egobp_dif_dn)


#tiff('ost_vic_dif_dn_bp.tiff', units="in", width=7, height=7, res=300, compression = 'lzw')
ggplot(egobp_dif_dn, aes(x= Dataset, y=Description, size=Count, color=p.adjust)) + 
  geom_point(alpha = 0.8) +
  scale_color_viridis() +
  theme_bw()
dev.off()


#Figure to paper
#upreg
egobp_ProtW_up$Dataset <- "3. Prw_contr"
egobp_ProtRel_up$Dataset <- "5. PrR_contr"
egobp_ProtW_dif_up$Dataset <- "4. Prw_dif"
egobp_ProtRel_dif_up$Dataset <- "6. PrR_dif"

egobp_ProtW_dif_up <- arrange(egobp_ProtW_dif_up, egobp_ProtW_dif_up$Count)
egobp_ProtW_dif_up1 <- egobp_ProtW_dif_up[c(1:20),]

egobp_up_all <- rbind(egobp_ProtW_up, egobp_ProtRel_up, egobp_ProtW_dif_up1, egobp_ProtRel_dif_up)
head(egobp_up_all)
egobp_up_all <- arrange(egobp_up_all, egobp_up_all$Count)

#tiff('ost_vic_full_up.tiff', units="in", width=8, height=7, res=300, compression = 'lzw')
ggplot(egobp_up_all, aes(x= Dataset, y=reorder(Description, Count), size=Count, color=p.adjust)) + 
  geom_point(alpha = 0.8) +
  scale_color_viridis() +
  theme_bw()
dev.off()

#down
egobp_ProtW_dn$Dataset <- "1. ProtW_con_d"
egobp_ProtRel_dn$Dataset <- "2. ProtR_con_d"
egobp_ProtRel_dif_dn$Dataset <- "3. ProtR_dif_d"


egobp_ProtRel_dn <- arrange(egobp_ProtRel_dn, egobp_ProtRel_dn$Count)
egobp_ProtRel_dn1 <- egobp_ProtRel_dn[c(1:20),]


egobp_dn_all <- rbind(egobp_ProtW_dn, egobp_ProtRel_dn1, egobp_ProtRel_dif_dn)
head(egobp_dn_all)

#tiff('ost_vic_full_dn.tiff', units="in", width=8, height=4, res=300, compression = 'lzw')
ggplot(egobp_dn_all, aes(x= Dataset, y=reorder(Description, Count), size=Count, color=p.adjust)) + 
  geom_point(alpha = 0.8) +
  scale_color_viridis() +
  theme_bw()
dev.off()


#Venn

unique_PrW_up <- VennDiagram::get.venn.partitions(list(contr=ProtW_up, diff=ProtW_dif_up))
unique_PrR_up <- VennDiagram::get.venn.partitions(list(contr=ProtRel_up, diff=ProtRel_dif_up))
unique_PrW_dn <- VennDiagram::get.venn.partitions(list(contr=ProtW_dn, diff=ProtW_dif_dn))
unique_PrR_dn <- VennDiagram::get.venn.partitions(list(contr=ProtRel_dn, diff=ProtRel_dif_dn))

venn.diagram(
  x = list(ProtW_up, ProtW_dif_up),
  category.names = c("Contr" , "Dif"),
  filename = '#unique_PrW_up.png',
  resolution = 300,
  print.mode=c("raw","percent"),
  output=F
)
venn.diagram(
  x = list(ProtRel_up, ProtRel_dif_up),
  category.names = c("Contr" , "Dif"),
  filename = '#unique_PrR_up.png',
  resolution = 300,
  print.mode=c("raw","percent"),
  output=F
)
venn.diagram(
  x = list(ProtW_dn, ProtW_dif_dn),
  category.names = c("Contr" , "Dif"),
  filename = '#unique_PrW_dn.png',
  resolution = 300,
  print.mode=c("raw","percent"),
  output=F
)
venn.diagram(
  x = list(ProtRel_dn, ProtRel_dif_dn),
  category.names = c("Contr" , "Dif"),
  filename = '#unique_PrR_dn.png',
  resolution = 300,
  print.mode=c("raw","percent"),
  output=F
)
