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
Tumor<-filt[,"C110rf95RELA">0]

### Scale data ###
ScaleData(Tumor)->Tumor

### Run MDS on Oncogene cluster ###
Tumor<-RunPCA(Tumor,dims=1:20)
Tumor<-RunTSNE(Tumor, dims=1:20)
Tumor<-RunUMAP(Tumor, dims=1:20)

### Visualize with MDS plots ###
pdf('Tumor_TSNE.pdf')
DimPlot(Tumor, reduction = "tsne", label = TRUE, pt.size = 0.75, label.size = 8)
pdf('Tumor_UMAP.pdf')
DimPlot(Tumor, reduction = "umap", label = TRUE, pt.size = 0.75, label.size = 8)
pdf('Tumor_PCA.pdf')
DimPlot(Tumor, reduction = "pca", label = TRUE, pt.size = 0.75, label.size = 8)

### Determine how many clusters should be used to group the cells ###
pdf('Tumor_Elbow.pdf')
ElbowPlot(Tumor, ndims=30)
pdf('Tumor_30Cluster_Heatmap.pdf')
DimHeatmap(Tumor, dims=1:30, cells=3000, balanced =TRUE)
Tumor <- JackStraw(Tumor, num.replicate=100)
pdf('Tumor_JackStraw.pdf')
Tumor <- ScoreJackStraw(Tumor, dims=1:20)
JackStrawPlot(Tumor, dims=1:20)

### Cluster the cells ###
Tumor <- FindNeighbors(Tumor, dims = 1:30)
Tumor <- FindClusters(Tumor, resolution = 0.35)

### Visualize the MDS plots ###
pdf('Tumor_TSNE.pdf')
DimPlot(Tumor, reduction = "tsne", label = TRUE, pt.size = 1.25, label.size = 8)
pdf('Tumor_UMAP.pdf')
DimPlot(Tumor, reduction = "umap", label = TRUE, pt.size = 1.25, label.size = 8)
pdf('Tumor_PCA.pdf')
DimPlot(Tumor, reduction = "pca", label = TRUE, pt.size = 1.25, label.size = 8)

### Identify Oncogene (C110rf95-RELA) cluster along with other ST-EPN signature genes ###
pdf('Tumor_C110rf95RELA_VlnPlot.pdf')
VlnPlot(Tumor, features = c("C110rf95RELA"), pt.size = 0)
pdf('Tumor_Lhx2_VlnPlot.pdf')
VlnPlot(Tumor, features = c("Lhx2"), pt.size = 0)
pdf('Tumor_Piezo1_VlnPlot.pdf')
VlnPlot(Tumor, features = c("Piezo1"), pt.size = 0)
pdf('Tumor_FeaturePlot.pdf')
FeaturePlot(Tumor, features= c("C110rf95RELA","Lhx2","Piezo1"))

### Identify marker genes and append info about gene descriptions ###
Tumor_markers <- FindAllMarkers(Tumor, only.pos = TRUE)
Tumor_mart <- useDataset("mmusculus_gene_ensembl", useMart(biomart = "ENSEMBL_MART_ENSEMBL"))
Tumor_genes <- Tumor_markers$gene
Tumor_descriptions <- getBM(attributes = c("external_gene_name", "description"), values = Tumor_genes, mart=Tumor_mart)
Tumor_out <- merge(as.data.frame(Tumor_markers), Tumor_descriptions, by.x= "gene", by.y = "external_gene_name")
write.csv(Tumor_out, file = "./Tumor_markers.csv")

### Identify the top 20 genes that distinguish a given cluster ###
Tumor_top20 <- Tumor_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
Tumor_genes <- Tumor_top20$gene
Tumor_top20descriptions <- getBM(attributes = c("external_gene_name", "description"),values=Tumor_genes,mart=Tumor_mart)
Tumor_top20out <- merge(as.data.frame(Tumor_top20), Tumor_descriptions, by.x= "gene", by.y="external_gene_name",all.x=TRUE)
write.csv(Tumor_top20out, file = "./Tumor_top20_markers.csv")

### Make a heatmap with the top marker genes from each cluster ###
pdf('Tumor_top20genes_heatmap.pdf')
DoHeatmap(Tumor, features = Tumor_top20$gene)

### Save Seurat Objects ###
saveRDS(Tumor,file="./Combined/Tumor_Cells/Combined_Tumor_scRNASeq_ClusteredSeurat.rds")