### Load the necessary packages for analysis###
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(SeuratData)
library(biomaRt)

workspace<-"/home/CAM/jhe/scRNASeq/Analysis/"
setwd(workspace)

filt <- readRDS("./Combined/Combined_FUS1_scRNASeq_Seurat.rds")

### Isolate C110rf95RELA-expressing clusters ###
Onc<-filt[,"cluster"=0]

### Scale data ###
ScaleData(Onc)->Onc

### Run MDS on oncogene cluster ###
Onc<-RunPCA(Onc,dims=1:20)
Onc<-RunTSNE(Onc, dims=1:20)
Onc<-RunUMAP(Onc, dims=1:20)

### Visualize with MDS plots ###
pdf('Onc_TSNE.pdf')
DimPlot(Onc, reduction = "tsne", label = TRUE, pt.size = 0.75, label.size = 8)
pdf('Onc_UMAP.pdf')
DimPlot(Onc, reduction = "umap", label = TRUE, pt.size = 0.75, label.size = 8)
pdf('Onc_PCA.pdf')
DimPlot(Onc, reduction = "pca", label = TRUE, pt.size = 0.75, label.size = 8)

### Determine how many clusters should be used to group the cells ###
pdf('Onc_Elbow.pdf')
ElbowPlot(Onc, ndims=30)
pdf('Onc_30Cluster_Heatmap.pdf')
DimHeatmap(Onc, dims=1:30, cells=3000, balanced =TRUE)
Onc <- JackStraw(Onc, num.replicate=100)
pdf('Onc_JackStraw.pdf')
Onc <- ScoreJackStraw(Onc, dims=1:20)
JackStrawPlot(Onc, dims=1:20)

### Cluster the cells ###
Onc <- FindNeighbors(Onc, dims = 1:30)
Onc <- FindClusters(Onc, resolution = 0.35)

### Visualize the MDS plots ###
pdf('onc_TSNE.pdf')
DimPlot(Onc, reduction = "tsne", label = TRUE, pt.size = 0.75, label.size = 8)
pdf('onc_UMAP.pdf')
DimPlot(Onc, reduction = "umap", label = TRUE, pt.size = 0.75, label.size = 8)
pdf('onc_PCA.pdf')
DimPlot(Onc, reduction = "pca", label = TRUE, pt.size = 0.75, label.size = 8)

### Identify oncogene (C110rf95-RELA) cluster along with other ST-EPN signature genes ###
pdf('Onc_C110rf95RELA_VlnPlot.pdf')
VlnPlot(Onc, features = c("C110rf95RELA"), pt.size = 0)
pdf('Onc_Lhx2_VlnPlot.pdf')
VlnPlot(Onc, features = c("Lhx2"), pt.size = 0)
pdf('Onc_Piezo1_VlnPlot.pdf')
VlnPlot(Onc, features = c("Piezo1"), pt.size = 0)
pdf('Onc_FeaturePlot.pdf')
FeaturePlot(Onc, features= c("C110rf95RELA","Lhx2","Piezo1"))

### Identify marker genes and append info about gene descriptions ###
Onc_markers <- FindAllMarkers(Onc, only.pos = TRUE)
Onc_mart <- useDataset("mmusculus_gene_ensembl", useMart(biomart = "ENSEMBL_MART_ENSEMBL"))
Onc_genes <- Onc_markers$gene
Onc_descriptions <- getBM(attributes = c("external_gene_name", "description"), values = Onc_genes, mart=Onc_mart)
Onc_out <- merge(as.data.frame(Onc_markers), Onc_descriptions, by.x= "gene", by.y = "external_gene_name")
write.csv(Onc_out, file = "./onc_markers.csv")

### Identify the top 20 genes that distinguish a given cluster ###
Onc_top20 <- Onc_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
Onc_genes <- Onc_top20$gene
Onc_top20descriptions <- getBM(attributes = c("external_gene_name", "description"),values=Onc_genes,mart=Onc_mart)
Onc_top20out <- merge(as.data.frame(Onc_top20), Onc_descriptions, by.x= "gene", by.y="external_gene_name",all.x=TRUE)
write.csv(Onc_top20out, file = "./Onc_top20_markers.csv")

### Make a heatmap with the top marker genes from each cluster ###
pdf('onc_top20genes_heatmap.pdf')
DoHeatmap(Onc, features = Onc_top20$gene)

### Save Seurat Objects ###
saveRDS(Onc,file="./Combined/Oncogene_Cluster/Combined_Oncogene_scRNASeq_ClusteredSeurat.rds")