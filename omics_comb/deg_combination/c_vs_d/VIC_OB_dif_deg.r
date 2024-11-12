#Osteoblasts
dda_ost <- data.frame(read.csv("dif_dda_ost_contVSost_diff.csv", header = TRUE)) 
dia_ost <- data.frame(read.csv("dif_dia_ost_contVSost_diff.csv", header = TRUE)) 
diaml_ost <- data.frame(read.csv("dif_diaml_ost_contVSost_diff.csv", header = TRUE)) 

dda_ost_up <- dda_ost[dda_ost$logFC > 1 & dda_ost$adj.P.Val < 0.05, 1] 
dda_ost_dn <- dda_ost[dda_ost$logFC < -1 & dda_ost$adj.P.Val < 0.05, 1] 

dia_ost_up <- dia_ost[dia_ost$logFC > 1 & dia_ost$adj.P.Val < 0.05, 1] 
dia_ost_dn <- dia_ost[dia_ost$logFC < -1 & dia_ost$adj.P.Val < 0.05, 1] 

diaml_ost_up <- diaml_ost[diaml_ost$logFC > 1 & diaml_ost$adj.P.Val < 0.05, 1] 
diaml_ost_dn <- diaml_ost[diaml_ost$logFC < -1 & diaml_ost$adj.P.Val < 0.05, 1] 




#Enrichment plot for upregulated
ProtW_ost_up <- c(dia_ost_up, diaml_ost_up, dda_ost_up)
ProtW_ost_up <- unique(ProtW_ost_up)
ProtRel_ost_up <- VennDiagram::get.venn.partitions(list(dia=dia_ost_up, diaml=diaml_ost_up, dda=dda_ost_up))
ProtRel_ost_up <- ProtRel_ost_up[1,5]
ProtRel_ost_up <- as.data.frame(ProtRel_ost_up)
ProtRel_ost_up <- ProtRel_ost_up$X1

library("org.Hs.eg.db")
library(rWikiPathways)
library(clusterProfiler)

ProtW_ost_up <- clusterProfiler::bitr(ProtW_ost_up,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
ProtRel_ost_up <- clusterProfiler::bitr(ProtRel_ost_up,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


egobp_ProtW_ost_up <- clusterProfiler::enrichGO(
  gene     = ProtW_ost_up[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05, 
  readable = TRUE)

egobp_ProtRel_ost_up <- clusterProfiler::enrichGO(
  gene     = ProtRel_ost_up[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE)


dotplot(egobp_ProtW_ost_up, showCategory = 20)
dotplot(egobp_ProtRel_ost_up, showCategory = 20)


library(viridis)
library(hrbrthemes)
library(ggplot2)

#Enrichment plot for downregulated
ProtW_ost_dn <- c(dia_ost_dn, diaml_ost_dn, dda_ost_dn)
ProtW_ost_dn <- unique(ProtW_ost_dn)
ProtRel_ost_dn <- VennDiagram::get.venn.partitions(list(dia=dia_ost_dn, diaml=diaml_ost_dn, dda=dda_ost_dn))
ProtRel_ost_dn <- ProtRel_ost_dn[1,5]
ProtRel_ost_dn <- as.data.frame(ProtRel_ost_dn)
ProtRel_ost_dn <- ProtRel_ost_dn$X1

ProtW_ost_dn <- clusterProfiler::bitr(ProtW_ost_dn,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
ProtRel_ost_dn <- clusterProfiler::bitr(ProtRel_ost_dn,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


egobp_ProtW_ost_dn <- clusterProfiler::enrichGO(
  gene     = ProtW_ost_dn[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05, 
  readable = TRUE)

egobp_ProtRel_ost_dn <- clusterProfiler::enrichGO(
  gene     = ProtRel_ost_dn[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE)

dotplot(egobp_ProtW_ost_dn, showCategory = 20)
dotplot(egobp_ProtRel_ost_dn, showCategory = 20)

egobp_ProtW_ost_dn <- data.frame(egobp_ProtW_ost_dn)
egobp_ProtRel_ost_dn <- data.frame(egobp_ProtRel_ost_dn)
#head(egobp_ProtW_ost_dn)



#########VICs

dda_VIC <- data.frame(read.csv("dda_vic_contVSvic_diff.csv", header = TRUE)) 
dia_VIC <- data.frame(read.csv("dia_vic_contVSvic_diff.csv", header = TRUE)) 
diaml_VIC <- data.frame(read.csv("diaml_vic_contVSvic_diff.csv", header = TRUE)) 
dda_VIC_up <- dda_VIC[dda_VIC$logFC > 1 & dda_VIC$adj.P.Val < 0.05, 1] 
dda_VIC_dn <- dda_VIC[dda_VIC$logFC < -1 & dda_VIC$adj.P.Val < 0.05, 1] 

dia_VIC_up <- dia_VIC[dia_VIC$logFC > 1 & dia_VIC$adj.P.Val < 0.05, 1] 
dia_VIC_dn <- dia_VIC[dia_VIC$logFC < -1 & dia_VIC$adj.P.Val < 0.05, 1] 

diaml_VIC_up <- diaml_VIC[diaml_VIC$logFC > 1 & diaml_VIC$adj.P.Val < 0.05, 1] 
diaml_VIC_dn <- diaml_VIC[diaml_VIC$logFC < -1 & diaml_VIC$adj.P.Val < 0.05, 1] 


#Enrichment plot for upregulated
ProtW_VIC_up <- c(dia_VIC_up, diaml_VIC_up, dda_VIC_up)
ProtW_VIC_up <- unique(ProtW_VIC_up)
ProtRel_VIC_up <- VennDiagram::get.venn.partitions(list(dia=dia_VIC_up, diaml=diaml_VIC_up, dda=dda_VIC_up))
ProtRel_VIC_up <- ProtRel_VIC_up[1,5]
ProtRel_VIC_up <- as.data.frame(ProtRel_VIC_up)
ProtRel_VIC_up <- ProtRel_VIC_up$X1

ProtW_VIC_up <- clusterProfiler::bitr(ProtW_VIC_up,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
ProtRel_VIC_up <- clusterProfiler::bitr(ProtRel_VIC_up,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


egobp_ProtW_VIC_up <- clusterProfiler::enrichGO(
  gene     = ProtW_VIC_up[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05, 
  readable = TRUE)

egobp_ProtRel_VIC_up <- clusterProfiler::enrichGO(
  gene     = ProtRel_VIC_up[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE)


dotplot(egobp_ProtW_VIC_up, showCategory = 20)
dotplot(egobp_ProtRel_VIC_up, showCategory = 20)

egobp_ProtW_VIC_up <- data.frame(egobp_ProtW_VIC_up)
egobp_ProtRel_VIC_up <- data.frame(egobp_ProtRel_VIC_up)

egobp_ProtW_VIC_up <- egobp_ProtW_VIC_up[,c(2,6,9)]
egobp_ProtW_VIC_up$Dataset <- "Prot. w.VIC"
egobp_ProtRel_VIC_up <- egobp_ProtRel_VIC_up[,c(2,6,9)]
egobp_ProtRel_VIC_up$Dataset <- "Prot. rel.VIC"

egobp_VIC_up <- rbind(egobp_ProtRel_VIC_up, egobp_ProtW_VIC_up)
head(egobp_VIC_up)

#tiff('vic_contrVSdif_up_bp.tiff', units="in", width=10, height=12, res=300, compression = 'lzw')
ggplot(egobp_VIC_up, aes(x= Dataset, y=Description, size=Count, color=p.adjust)) + 
  geom_point(alpha = 0.8) +
  scale_color_viridis() +
  theme_bw()
dev.off()

#Enrichment plot for downregulated
ProtW_VIC_dn <- c(dia_VIC_dn, diaml_VIC_dn, dda_VIC_dn)
ProtW_VIC_dn <- unique(ProtW_VIC_dn)
ProtRel_VIC_dn <- VennDiagram::get.venn.partitions(list(dia=dia_VIC_dn, diaml=diaml_VIC_dn, dda=dda_VIC_dn))
ProtRel_VIC_dn <- ProtRel_VIC_dn[1,5]
ProtRel_VIC_dn <- as.data.frame(ProtRel_VIC_dn)
ProtRel_VIC_dn <- ProtRel_VIC_dn$X1

ProtW_VIC_dn <- clusterProfiler::bitr(ProtW_VIC_dn,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
ProtRel_VIC_dn <- clusterProfiler::bitr(ProtRel_VIC_dn,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


egobp_ProtW_VIC_dn <- clusterProfiler::enrichGO(
  gene     = ProtW_VIC_dn[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05, 
  readable = TRUE)

egobp_ProtRel_VIC_dn <- clusterProfiler::enrichGO(
  gene     = ProtRel_VIC_dn[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE)


egobp_ProtW_VIC_dn <- data.frame(egobp_ProtW_VIC_dn)
egobp_ProtRel_VIC_dn <- data.frame(egobp_ProtRel_VIC_dn)

egobp_ProtRel_VIC_dn <- egobp_ProtRel_VIC_dn[,c(2,6,9)]
egobp_ProtRel_VIC_dn$Dataset <- "ProtRel_VIC_dn"
egobp_ProtW_VIC_dn <- egobp_ProtW_VIC_dn[,c(2,6,9)]
egobp_ProtW_VIC_dn$Dataset <- "ProtW_VIC_dn"

egobp_VIC_dn <- rbind(egobp_ProtRel_VIC_dn, egobp_ProtW_VIC_dn)
head(egobp_VIC_dn)


#tiff('vic_contrVSdif_dn_bp.tiff', units="in", width=7, height=7, res=300, compression = 'lzw')
ggplot(egobp_VIC_dn, aes(x= Dataset, y=Description, size=Count, color=p.adjust)) + 
  geom_point(alpha = 0.8) +
  scale_color_viridis() +
  theme_bw()
dev.off()