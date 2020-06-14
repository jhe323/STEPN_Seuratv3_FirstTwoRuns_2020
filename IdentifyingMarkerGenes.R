### Load the necessary packages for analysis###
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(SeuratData)
library(biomaRt)

run1<-"/home/CAM/jhe/scRNASeq/Analysis/Run1"
setwd(run1)
r1 <- readRDS("./Run1/Run1_FUS1_scRNASeq_Seurat.rds")

run2<-"/home/CAM/jhe/scRNASeq/Analysis/Run2"
setwd(run2)
r2 <- readRDS("./Run2/Run2_FUS1_scRNASeq_Seurat.rds")

workspace<-"/home/CAM/jhe/scRNASeq/Analysis/"
setwd(workspace)

filt <- readRDS("./Combined/Combined_FUS1_scRNASeq_Seurat.rds")

### Identify marker genes and append info about gene descriptions ###
markers <- FindAllMarkers(filt, only.pos = TRUE)
mart <- useDataset("mmusculus_gene_ensembl", useMart(biomart = "ENSEMBL_MART_ENSEMBL"))
genes <- markers$gene
descriptions <- getBM(attributes = c("external_gene_name", "description"), values = genes, mart = mart)
out <- merge(as.data.frame(markers), descriptions, by.x= "gene", by.y = "external_gene_name")
write.csv(out, file = "./markers.csv")

r1_markers <- FindAllMarkers(r1, only.pos = TRUE)
r1_mart <- useDataset("mmusculus_gene_ensembl", useMart(biomart = "ENSEMBL_MART_ENSEMBL"))
r1_genes <- r1_markers$gene
r1_descriptions <- getBM(attributes = c("external_gene_name", "description"), values = r1_genes, mart=r1_mart)
r1_out <- merge(as.data.frame(r1_markers), r1_descriptions, by.x= "gene", by.y = "external_gene_name")
write.csv(r1_out, file = "./r1_markers.csv")

r2_markers <- FindAllMarkers(r2, only.pos = TRUE)
r2_mart <- useDataset("mmusculus_gene_ensembl", useMart(biomart = "ENSEMBL_MART_ENSEMBL"))
r2_genes <- r2_markers$gene
r2_descriptions <- getBM(attributes = c("external_gene_name", "description"), values = r2_genes, mart=r2_mart)
r2_out <- merge(as.data.frame(r2_markers), r2_descriptions, by.x= "gene", by.y = "external_gene_name")
write.csv(r2_out, file = "./r2_markers.csv")


### Identify the top 20 genes that distinguish a given cluster ###
r1_top20 <- r1_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
r1_genes <- r1_top20$gene
r1_top20descriptions <- getBM(attributes = c("external_gene_name","description"), values=r1_genes, mart=r1_mart)
r1_top20out <- merge(as.data.frame(r1_top20), r1_descriptions, by.x= "gene", by.y="external_gene_name",all.x=TRUE)
write.csv(r1_top20out, file = "./r1_top20_markers.csv")

r2_top20 <- r2_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
r2_genes <- r2_top20$gene
r2_top20descriptions <- getBM(attributes = c("external_gene_name", "description"),values=r2_genes,mart=r2_mart)
r2_top20out <- merge(as.data.frame(r2_top20), r2_descriptions, by.x= "gene", by.y = "external_gene_name",all.x=TRUE)
write.csv(r2_top20out, file = "./r2_top20_markers.csv")

top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
genes <- top20$gene
top20descriptions <- getBM(attributes = c("external_gene_name", "description"), values = genes, mart = mart)
top20out <- merge(as.data.frame(top20), descriptions, by.x= "gene", by.y = "external_gene_name",all.x = TRUE)
write.csv(top20out, file = "./top20_markers.csv")


### Make a heatmap with the top marker genes from each cluster ###
pdf('Combined_top20genes_heatmap.pdf')
DoHeatmap(filt, features = top20$gene)

pdf('r1_top20genes_heatmap.pdf')
DoHeatmap(r1, features = r2_top20$gene)

pdf('r2_top20genes_heatmap.pdf')
DoHeatmap(r2, features = r2_top20$gene)
