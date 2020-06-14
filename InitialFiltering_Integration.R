### Load the necessary packages for analysis###
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(SeuratData)
library(biomaRt)

### Load the cell cycle genes dataset from Seurat for cell cycle sorting later on###
data("cc.genes")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

### Set the proper working directories containing the filtered_gene_bc_matrices for analysis###
Run1<-"/home/CAM/jhe/scRNASeq/Run1/"
Run2<-"/home/CAM/jhe/scRNASeq/Run2/"
workspace<-"/home/CAM/jhe/scRNASeq/Analysis/"

### Load the respective files you need and create Seurat objects for each###
setwd(Run1)
Read10X(data.dir="./filtered_feature_bc_matrix/")->r1
r1<-CreateSeuratObject(counts=r1,project="STEPN_Run1",min.cells=3,min.features=200)
setwd(Run2)
Read10X(data.dir="./filtered_feature_bc_matrix/")->r2
r2<-CreateSeuratObject(counts=r2,project="STEPN_Run2",min.cells=3,min.features=200)

### Note: Depending on how many cells you sequenced, these objects may get pretty large (>1 gb) and slow down computation###

setwd(workspace)

### Add mitochondrial content and designate a tag for each run into metadata ###
r1[["percent.mt"]] <- PercentageFeatureSet(r1, pattern = "^mt-")
r1[["replicate"]] <- "1"
r2[["percent.mt"]] <- PercentageFeatureSet(r2, pattern = "^mt-")
r2[["replicate"]] <- "2"

### Normalize the data ###
r1 <- NormalizeData(r1)
r2 <- NormalizeData(r2)

### Find variable features (DNA expression) unique to each run of cells ###
r1 <- FindVariableFeatures(r1, selection.method = "vst", nfeatures=2000)
r2 <- FindVariableFeatures(r2, selection.method = "vst", nfeatures=2000)

### Label each cell for their respective cell cycle state ###
r1 <- CellCycleScoring(r1, s.features = s.genes, g2m.features = g2m.genes)
r2 <- CellCycleScoring(r2, s.features = s.genes, g2m.features = g2m.genes)

### Find integration anchors in the datasets and integrate the two runs together ###
df_anchors <- FindIntegrationAnchors(object.list=c(r1,r2), dims=1:30)
df <- IntegrateData(anchorset= df_anchors, dims=1:30)
df <- ScaleData(df)

ScaleData(r1)->r1
ScaleData(r2)->r2

### Make a Violin Plot of pre-filtered metadata ###
pdf('Run1_RawQC.pdf')
VlnPlot(r1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pdf('Run2_RawQC.pdf')
VlnPlot(r2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pdf('Combined_Filtered_QC.pdf')
VlnPlot(df, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

### Run multi-dimensional scaling (PCA, tSNE, UMAP) of the combined dataset with the unfiltered data to see how it looks ###
df <- RunPCA(df, dims = 1:20)
df <- RunTSNE(df, dims =1:20)
df<- RunUMAP(df, dims=1:20)

r1<-RunPCA(r1,dims=1:20)
r1<-RunTSNE(r1, dims=1:20)
r1<-RunUMAP(r1, dims=1:20)

r2<-RunPCA(r2,dims=1:20)
r2<-RunTSNE(r2, dims=1:20)
r2<-RunUMAP(r2, dims=1:20)

### Check for batch effects between the two runs ###
pdf('Combined_Batch.df')
DimPlot(df, reduction = 'umap', group.by = 'replicate')

### Make PCA, tSNE, and UMAP plots (unfiltered) ###
DimPlot(df, reduction = "tsne", label = TRUE, pt.size = 0.75, label.size = 8)
DimPlot(df, reduction = "umap", label = TRUE, pt.size = 0.75, label.size = 8)
DimPlot(df, reduction = "pca", label = TRUE, pt.size = 0.75, label.size = 8)

### Filter out cells with high mitochondrial RNA content based on how MDS looks ###
filt <- subset(df, subset = percent.mt <= 20 & nFeature_RNA <= 6000 & nCount_RNA <= 25000)
r1 <- subset(r1, subset = percent.mt <= 20 & nFeature_RNA <= 6000 & nCount_RNA <= 25000)
r2 <- subset(r2, subset = percent.mt <= 20 & nFeature_RNA <= 6000 & nCount_RNA <= 25000)

### Visualize filtered data with a violin plot ###
pdf('Run1_FilteredQC.pdf')
VlnPlot(r1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pdf('Run2_FilteredQC.pdf')
VlnPlot(r2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pdf('Combined_FilteredQC.pdf')
VlnPlot(filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

### Save this Seurat Object as a checkpoint ###
saveRDS(filt, file = "./Combined/Combined_FUS1_scRNASeq_Seurat.rds")
saveRDS(r1, file= "./Run1/Run1_FUS1_scRNASeq_Seurat.rds")
saveRDS(r2,file="./Run2/Run2_FUS1_scRNASeq_Seurat.rds")