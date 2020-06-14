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

### Determine how many clusters should be used to group the cells ###
pdf('Combined_Elbow.pdf')
ElbowPlot(filt, ndims=30)
pdf('Combined_30Cluster_Heatmap.pdf')
DimHeatmap(filt, dims=1:30, cells=3000, balanced =TRUE)
filt <- JackStraw(filt, num.replicate=100)
pdf('Combined_JackStraw.pdf')
filt <- ScoreJackStraw(filt, dims=1:20)
JackStrawPlot(filt, dims=1:20)

### Cluster the cells ###
filt <- FindNeighbors(filt, dims = 1:20)
filt <- FindClusters(filt, resolution = 0.35)

r1 <- FindNeighbors(r1, dims = 1:30)
r1 <- FindClusters(r1, resolution = 0.35)

r2 <- FindNeighbors(r2, dims = 1:30)
r2 <- FindClusters(r2, resolution = 0.35)

### Visualize with MDS plots ###
pdf('Combined_TSNE.pdf')
DimPlot(filt, reduction = "tsne", label = TRUE, pt.size = 0.75, label.size = 8)
pdf('Combined_UMAP.pdf')
DimPlot(filt, reduction = "umap", label = TRUE, pt.size = 0.75, label.size = 8)
pdf('Combined_PCA.pdf')
DimPlot(filt, reduction = "pca", label = TRUE, pt.size = 0.75, label.size = 8)

pdf('r1_TSNE.pdf')
DimPlot(r1, reduction = "tsne", label = TRUE, pt.size = 0.75, label.size = 8)
pdf('r1_UMAP.pdf')
DimPlot(r1, reduction = "umap", label = TRUE, pt.size = 0.75, label.size = 8)
pdf('r1_PCA.pdf')
DimPlot(r1, reduction = "pca", label = TRUE, pt.size = 0.75, label.size = 8)

pdf('r2_TSNE.pdf')
DimPlot(r2, reduction = "tsne", label = TRUE, pt.size = 0.75, label.size = 8)
pdf('r2_UMAP.pdf')
DimPlot(r2, reduction = "umap", label = TRUE, pt.size = 0.75, label.size = 8)
pdf('r2_PCA.pdf')
DimPlot(r2, reduction = "pca", label = TRUE, pt.size = 0.75, label.size = 8)

### Identify oncogene (C110rf95-RELA) cluster along with other ST-EPN signature genes (Lhx2, Piezo1, Rela) ###
pdf('Combined_C110rf95RELA_VlnPlot.pdf')
VlnPlot(filt, features = c("C110rf95RELA"), pt.size = 0)
pdf('Combined_Lhx2_VlnPlot.pdf')
VlnPlot(filt, features = c("Lhx2"), pt.size = 0)
pdf('Combined_piezo1_VlnPlot.pdf')
VlnPlot(filt, features = c("Piezo1"), pt.size = 0)
FeaturePlot(filt, features= c("C110rf95RELA","Lhx2","Piezo1"),pt.size=0.50)

pdf('r1_C110rf95RELA_VlnPlot.pdf')
VlnPlot(r1, features = c("C110rf95RELA"), pt.size = 0)
pdf('r1_Lhx2_VlnPlot.pdf')
VlnPlot(r1, features = c("Lhx2"), pt.size = 0)
pdf('r1_Piezo1_VlnPlot.pdf')
VlnPlot(r1, features = c("Piezo1"), pt.size = 0)
pdf('r1_FeaturePlot.pdf')
FeaturePlot(r1, features= c("C110rf95RELA","Lhx2","Piezo1"),pt.size=0.50)

pdf('r2_C110rf95RELA_VlnPlot.pdf')
VlnPlot(r2, features = c("C110rf95RELA"), pt.size = 0)
pdf('r2_Lhx2_VlnPlot.pdf')
VlnPlot(r2, features = c("Lhx2"), pt.size = 0)
pdf('r2_Piezo1_VlnPlot.pdf')
VlnPlot(r2, features = c("Piezo1"), pt.size = 0)
pdf('r2_FeaturePlot.pdf')
FeaturePlot(r2, features= c("C110rf95RELA","Lhx2","Piezo1"),pt.size=0.50)

### Save this Seurat Object as a checkpoint ###
saveRDS(filt, file = "./Combined/Combined_FUS1_scRNASeq_Seurat.rds")
saveRDS(r1, file= "./Run1/Run1_FUS1_scRNASeq_Seurat.rds")
saveRDS(r2,file="./Run2/Run2_FUS1_scRNASeq_Seurat.rds")