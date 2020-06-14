### Load the necessary packages for analysis###
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(SeuratData)
library(biomaRt)

### Load the cell cycle genes dataset from Seurat for cell cycle sorting later on###
#data("cc.genes")
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes

### Set the proper working directories containing the filtered_gene_bc_matrices for analysis###
Run1<-"/home/CAM/jhe/scRNASeq/Run1/"
Run2<-"/home/CAM/jhe/scRNASeq/Run2/"
workspace<-"/home/CAM/jhe/scRNASeq/Analysis/"

### Load the respective files you need and create Seurat objects for each###
#setwd(Run1)
#Read10X(data.dir="./filtered_feature_bc_matrix/")->r1
#r1<-CreateSeuratObject(counts=r1,project="STEPN_Run1",min.cells=3,min.features=200)
#setwd(Run2)
#Read10X(data.dir="./filtered_feature_bc_matrix/")->r2
#r2<-CreateSeuratObject(counts=r2,project="STEPN_Run2",min.cells=3,min.features=200)

### Note: Depending on how many cells you sequenced, these objects may get pretty large (>1 gb) and slow down computation###

setwd(workspace)

### Add mitochondrial content and designate a tag for each run into metadata ###
#r1[["percent.mt"]] <- PercentageFeatureSet(r1, pattern = "^mt-")
#r1[["replicate"]] <- "1"
#r2[["percent.mt"]] <- PercentageFeatureSet(r2, pattern = "^mt-")
#r2[["replicate"]] <- "2"

### Normalize the data ###
#r1 <- NormalizeData(r1)
#r2 <- NormalizeData(r2)

### Find variable features (DNA expression) unique to each run of cells ###
#r1 <- FindVariableFeatures(r1, selection.method = "vst", nfeatures=2000)
#r2 <- FindVariableFeatures(r2, selection.method = "vst", nfeatures=2000)

### Label each cell for their respective cell cycle state ###
#r1 <- CellCycleScoring(r1, s.features = s.genes, g2m.features = g2m.genes)
#r2 <- CellCycleScoring(r2, s.features = s.genes, g2m.features = g2m.genes)

### Find integration anchors in the datasets and integrate the two runs together ###
#df_anchors <- FindIntegrationAnchors(object.list=c(r1,r2), dims=1:30)
#df <- IntegrateData(anchorset= df_anchors, dims=1:30)
#df <- ScaleData(df)

#filt <- readRDS("./Combined/Combined_FUS1_scRNASeq_Seurat.rds")
#r1 <- readRDS("./Run1/Run1_FUS1_scRNASeq_Seurat.rds")
#r2 <- readRDS("./Run2/Run2_FUS1_scRNASeq_Seurat.rds")

### Make a Violin Plot of pre-filtered metadata ###
#pdf('Run1_RawQC.pdf')
#VlnPlot(r1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#pdf('Run2_RawQC.pdf')
#VlnPlot(r2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(df, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

### Run multi-dimensional scaling (PCA, tSNE, UMAP) of the combined dataset with the unfiltered data to see how it looks ###
#df <- RunPCA(df, dims = 1:30)
#df <- RunTSNE(df, dims =1:30)
#df<- RunUMAP(df, dims=1:30)

### Make PCA, tSNE, and UMAP plots (unfiltered) ###
#DimPlot(df, reduction = "tsne", label = TRUE, pt.size = 1.25, label.size = 8)
#DimPlot(df, reduction = "umap", label = TRUE, pt.size = 1.25, label.size = 8)
#DimPlot(df, reduction = "pca", label = TRUE, pt.size = 1.25, label.size = 8)

### Filter out cells with high mitochondrial RNA content based on how MDS looks ###
#filt <- subset(df, subset = percent.mt <= 20 & nFeature_RNA <= 6000 & nCount_RNA <= 25000)
#r1 <- subset(r1, subset = percent.mt <= 20 & nFeature_RNA <= 6000 & nCount_RNA <= 25000)
#r2 <- subset(r2, subset = percent.mt <= 20 & nFeature_RNA <= 6000 & nCount_RNA <= 25000)

### Visualize filtered data with a violin plot ###
#pdf('Run1_FilteredQC.pdf')
#VlnPlot(r1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#pdf('Run2_FilteredQC.pdf')
#VlnPlot(r2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#pdf('Combined_FilteredQC.pdf')
#VlnPlot(filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

### Determine how many clusters should be used to group the cells ###
#pdf('Combined_Elbow.pdf')
#ElbowPlot(filt, ndims=30)
#pdf('Combined_30Cluster_Heatmap.pdf')
#DimHeatmap(filt, dims=1:30, cells=3000, balanced =TRUE)
##filt <- JackStraw(filt, num.replicate=100)
#pdf('Combined_JackStraw.pdf')
#filt <- ScoreJackStraw(filt, dims=1:20)
#JackStrawPlot(filt, dims=1:20)

### Cluster the cells ###
#filt <- FindNeighbors(filt, dims = 1:20)
#filt <- FindClusters(filt, resolution = 0.35)

### Visualize with MDS plots ###
#pdf('Combined_TSNE.pdf')
#DimPlot(filt, reduction = "tsne", label = TRUE, pt.size = 1.25, label.size = 8)
#pdf('Combined_UMAP.pdf')
#DimPlot(filt, reduction = "umap", label = TRUE, pt.size = 1.25, label.size = 8)
#pdf('Combined_PCA.pdf')
#DimPlot(filt, reduction = "pca", label = TRUE, pt.size = 1.25, label.size = 8)

### Identify oncogene (C110rf95-RELA) cluster along with other ST-EPN signature genes (Lhx2, Piezo1, Rela) ###
#pdf('Combined_C110rf95RELA_VlnPlot.pdf')
#VlnPlot(filt, features = c("C110rf95RELA"), pt.size = 0)
#pdf('Combined_Lhx2_VlnPlot.pdf')
#VlnPlot(filt, features = c("Lhx2"), pt.size = 0)
#pdf('Combined_piezo1_VlnPlot.pdf')
#VlnPlot(filt, features = c("Piezo1"), pt.size = 0)
##FeaturePlot(filt, features= c("C110rf95RELA","Lhx2","Piezo1"))

### Identify marker genes and append info about gene descriptions ###
#markers <- FindAllMarkers(filt, only.pos = TRUE)
#mart <- useDataset("mmusculus_gene_ensembl", useMart(biomart = "ENSEMBL_MART_ENSEMBL"))
#genes <- markers$gene
#descriptions <- getBM(attributes = c("external_gene_name", "description"), values = genes, mart = mart)
#out <- merge(as.data.frame(markers), descriptions, by.x= "gene", by.y = "external_gene_name")
#write.csv(out, file = "./markers.csv")

### Identify the top 20 genes that distinguish a given cluster ###
#top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
#genes <- top20$gene
#top20descriptions <- getBM(attributes = c("external_gene_name", "description"), values = genes, mart = mart)
#top20out <- merge(as.data.frame(top20), descriptions, by.x= "gene", by.y = "external_gene_name",all.x = TRUE)
#write.csv(top20out, file = "./top20_markers.csv")

### Make a heatmap with the top marker genes from each cluster ###
#pdf('Combined_top20genes_heatmap.pdf')
#DoHeatmap(filt, features = top20$gene)

### Save this Seurat Object as a checkpoint ###
#saveRDS(filt, file = "./Combined/Combined_FUS1_scRNASeq_Seurat.rds")
#saveRDS(r1, file= "./Run1/Run1_FUS1_scRNASeq_Seurat.rds")
#saveRDS(r2,file="./Run2/Run2_FUS1_scRNASeq_Seurat.rds")

### Use this function to re-open a saved RDS Seurat Object ###
### filt <- readRDS("./filtered_SeuratObject.rds")

### Isolate C110rf95RELA-expressing clusters ###
#Onc<-filt[,"cluster"=0]

### Scale data ###
#ScaleData(Onc)->Onc
#ScaleData(r1)->r1
#ScaleData(r2)->r2

### Run MDS on oncogene cluster ###
#Onc<-RunPCA(Onc,dims=1:20)
#Onc<-RunTSNE(Onc, dims=1:20)
#Onc<-RunUMAP(Onc, dims=1:20)

#r1<-RunPCA(r1,dims=1:20)
#r1<-RunTSNE(r1, dims=1:20)
#r1<-RunUMAP(r1, dims=1:20)

#r2<-RunPCA(r2,dims=1:20)
#r2<-RunTSNE(r2, dims=1:20)
#r2<-RunUMAP(r2, dims=1:20)

### Cluster the cells ###
#r1 <- FindNeighbors(r1, dims = 1:30)
#r1 <- FindClusters(r1, resolution = 0.35)

#r2 <- FindNeighbors(r2, dims = 1:30)
#r2 <- FindClusters(r2, resolution = 0.35)

#Onc <- FindNeighbors(Onc, dims = 1:30)
#Onc <- FindClusters(Onc, resolution = 0.35)

### Visualize the MDS plots ###
#pdf('r1_TSNE.pdf')
#DimPlot(r1, reduction = "tsne", label = TRUE, pt.size = 1.25, label.size = 8)
#pdf('r1_UMAP.pdf')
#DimPlot(r1, reduction = "umap", label = TRUE, pt.size = 1.25, label.size = 8)
#pdf('r1_PCA.pdf')
#DimPlot(r1, reduction = "pca", label = TRUE, pt.size = 1.25, label.size = 8)

#pdf('r2_TSNE.pdf')
#DimPlot(r2, reduction = "tsne", label = TRUE, pt.size = 1.25, label.size = 8)
#pdf('r2_UMAP.pdf')
#DimPlot(r2, reduction = "umap", label = TRUE, pt.size = 1.25, label.size = 8)
#pdf('r2_PCA.pdf')
#DimPlot(r2, reduction = "pca", label = TRUE, pt.size = 1.25, label.size = 8)

#pdf('onc_TSNE.pdf')
#DimPlot(Onc, reduction = "tsne", label = TRUE, pt.size = 1.25, label.size = 8)
#pdf('onc_UMAP.pdf')
#DimPlot(Onc, reduction = "umap", label = TRUE, pt.size = 1.25, label.size = 8)
#pdf('onc_PCA.pdf')
#DimPlot(Onc, reduction = "pca", label = TRUE, pt.size = 1.25, label.size = 8)

### Identify oncogene (C110rf95-RELA) cluster along with other ST-EPN signature genes ###
#pdf('r1_C110rf95RELA_VlnPlot.pdf')
#VlnPlot(r1, features = c("C110rf95RELA"), pt.size = 0)
#pdf('r1_Lhx2_VlnPlot.pdf')
#VlnPlot(r1, features = c("Lhx2"), pt.size = 0)
#pdf('r1_Piezo1_VlnPlot.pdf')
#VlnPlot(r1, features = c("Piezo1"), pt.size = 0)
#pdf('r1_FeaturePlot.pdf')
#FeaturePlot(r1, features= c("C110rf95RELA","Lhx2","Piezo1"))

#pdf('r2_C110rf95RELA_VlnPlot.pdf')
#VlnPlot(r2, features = c("C110rf95RELA"), pt.size = 0)
#pdf('r2_Lhx2_VlnPlot.pdf')
#VlnPlot(r2, features = c("Lhx2"), pt.size = 0)
#pdf('r2_Piezo1_VlnPlot.pdf')
#VlnPlot(r2, features = c("Piezo1"), pt.size = 0)
#pdf('r2_FeaturePlot.pdf')
#FeaturePlot(r2, features= c("C110rf95RELA","Lhx2","Piezo1"))

#pdf('Onc_C110rf95RELA_VlnPlot.pdf')
#VlnPlot(Onc, features = c("C110rf95RELA"), pt.size = 0)
#pdf('Onc_Lhx2_VlnPlot.pdf')
#VlnPlot(Onc, features = c("Lhx2"), pt.size = 0)
#pdf('Onc_Piezo1_VlnPlot.pdf')
#VlnPlot(Onc, features = c("Piezo1"), pt.size = 0)
#pdf('Onc_FeaturePlot.pdf')
#FeaturePlot(Onc, features= c("C110rf95RELA","Lhx2","Piezo1"))

### Identify marker genes and append info about gene descriptions ###
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

Onc_markers <- FindAllMarkers(Onc, only.pos = TRUE)
Onc_mart <- useDataset("mmusculus_gene_ensembl", useMart(biomart = "ENSEMBL_MART_ENSEMBL"))
Onc_genes <- Onc_markers$gene
Onc_descriptions <- getBM(attributes = c("external_gene_name", "description"), values = Onc_genes, mart=Onc_mart)
Onc_out <- merge(as.data.frame(Onc_markers), Onc_descriptions, by.x= "gene", by.y = "external_gene_name")
write.csv(Onc_out, file = "./onc_markers.csv")

### Identify the top 20 genes that distinguish a given cluster ###
r1_top20 <- r1_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
r1_genes <- r1_top20$gene
r1_top20descriptions <- getBM(attributes = c("external_gene_name","description"), values=r1_genes, mart=r1_mart)

r2_top20 <- r2_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
r2_genes <- r2_top20$gene
r2_top20descriptions <- getBM(attributes = c("external_gene_name", "description"),values=r2_genes,mart=r2_mart)

r1_top20out <- merge(as.data.frame(r1_top20), r1_descriptions, by.x= "gene", by.y="external_gene_name",all.x=TRUE)
write.csv(r1_top20out, file = "./r1_top20_markers.csv")

r2_top20out <- merge(as.data.frame(r2_top20), r2_descriptions, by.x= "gene", by.y = "external_gene_name",all.x=TRUE)
write.csv(r2_top20out, file = "./r2_top20_markers.csv")

Onc_top20 <- Onc_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
Onc_genes <- Onc_top20$gene
Onc_top20descriptions <- getBM(attributes = c("external_gene_name", "description"),values=Onc_genes,mart=Onc_mart)
Onc_top20out <- merge(as.data.frame(Onc_top20), Onc_descriptions, by.x= "gene", by.y="external_gene_name",all.x=TRUE)
write.csv(Onc_top20out, file = "./Onc_top20_markers.csv")

### Make a heatmap with the top marker genes from each cluster ###
pdf('r1_top20genes_heatmap.pdf')
DoHeatmap(r1, features = r2_top20$gene)

pdf('r2_top20genes_heatmap.pdf')
DoHeatmap(r2, features = r2_top20$gene)

pdf('onc_top20genes_heatmap.pdf')
DoHeatmap(Onc, features = Onc_top20$gene)

### Save Seurat Objects ###
saveRDS(Onc,file="./Combined_Oncogene_scRNASeq_ClusteredSeurat.rds")
#saveRDS(r1, file= "./Run1/Run1_FUS1_scRNASeq_Seurat.rds")
#saveRDS(r2,file="./Run2/Run2_FUS1_scRNASeq_Seurat.rds")

