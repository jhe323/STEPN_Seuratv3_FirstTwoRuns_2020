### Load the necessary packages for analysis###
onc*library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(SeuratData)
library(biomaRt)
library(zinbwave)
library(scRNAseq)
library(matrixStats)
library(magrittr)
library(SingleCellExperiment)
library(stringr)

### Load the cell cycle genes dataset from Seurat for cell cycle sorting later on###
data("cc.genes")
s.genes <- cc.genes$s.genes
s.genes <- str_to_title(s.genes)
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- str_to_title(g2m.genes)

### Set the proper working directories containing the filtered_gene_bc_matrices for analysis###
#Run1<-"/home/CAM/jhe/scRNASeq/Run1/"
#Run2<-"/home/CAM/jhe/scRNASeq/Run2/"
#workspace<-"/home/CAM/jhe/scRNASeq/Analysis/"

### Load the respective files you need and create Seurat objects for each###
#setwd(Run1)
Read10X(data.dir="../Run1/filtered_feature_bc_matrix/")->r1
r1<-CreateSeuratObject(counts=r1,project="STEPN_Run1",min.cells=3,min.features=200)
#setwd(Run2)
Read10X(data.dir="../Run2/filtered_feature_bc_matrix/Seurat_Files/")->r2
r2<-CreateSeuratObject(counts=r2,project="STEPN_Run2",min.cells=3,min.features=200)

### Note: Depending on how many cells you sequenced, these objects may get pretty large (>1 gb) and slow down computation###

#setwd(workspace)

### Add mitochondrial content and designate a tag for each run into metadata ###
r1[["percent.mt"]] <- PercentageFeatureSet(r1, pattern = "^mt-")
r1[["replicate"]] <- "1"
r2[["percent.mt"]] <- PercentageFeatureSet(r2, pattern = "^mt-")
r2[["replicate"]] <- "2"

### Label each cell for their respective cell cycle state ###
r1 <- CellCycleScoring(r1, s.features = s.genes, g2m.features = g2m.genes)
r2 <- CellCycleScoring(r2, s.features = s.genes, g2m.features = g2m.genes)

### Convert Seurat object to a SingleCellExperiment Object ###
SingleCellExperiment(r1[["RNA"]]@data)->r1.sce
SingleCellExperiment(r2[["RNA"]]@data)->r2.sce

### ZINB WaVE Gene Filtering of Runs 1 and 2 and set up the assay name for zinbwave normalization ###
filter.1 <- rowSums(assay(r1.sce)>5)>5
table(filter.1)
#filter.1
r1.sce <- r1.sce[filter.1,]
rowVars(as.matrix(assay(r1.sce))) -> vars.1
names(vars.1) <- rownames(r1.sce)
vars.1 <- sort(vars.1, decreasing = TRUE)
head(vars.1)
r1.sce <- r1.sce[names(vars.1)[1:100],]
assayNames(r1.sce)[1] <- "counts"

filter.2 <- rowSums(assay(r2.sce)>5)>5
table(filter.2)
#filter.2
r2.sce <- r2.sce[filter.2,]
rowVars(as.matrix(assay(r2.sce))) -> vars.2
names(vars.2) <- rownames(r2.sce)
vars.2 <- sort(vars.2, decreasing = TRUE)
head(vars.2)
r2.sce <- r2.sce[names(vars.2)[1:100],]
assayNames(r2.sce)[1] <- "counts"

### Convert assay to matrix and run zinbwave normalization on runs 1 and 2 ###
as.matrix(assay(r1.sce))->assay(r1.sce)
ifelse(assay(r1.sce)==".",0,assay(r1.sce))->assay(r1.sce)
as.matrix(assay(r2.sce))->assay(r2.sce)
ifelse(assay(r2.sce)==".",0,assay(r2.sce))->assay(r2.sce)

#r1_zinb <- zinbwave(r1.sce, K = 20, epsilon=1000000, normalizedValues=T, residuals=T, verbose = T)
#r2_zinb <- zinbwave(r2.sce, K = 20, epsilon=1000000, normalizedValues=T, residuals=T, verbose = T)

r1_surf <- zinbsurf(r1.sce, K=2, epsilon = 1000, verbose = T, prop_fit = 0.5)

r1_surf

r2_surf <- zinbsurf(r2.sce, K=2, epsilon = 1000, verbose = T, prop_fit = 0.5)

r2_surf

W <- reducedDim(r1_surf)
X <- reducedDim(r2_surf)

r1_zinb
r2_zinb
write.table(W.1,"Run1_ZinbWave_W.xlsx",quote=F,sep="\t",col.names=T,row.names=T)
write.table(W.2,"Run2_ZinbWave_W.xlsx",quote=F,sep="\t",col.names=T,row.names=T)

### Visualize the W Matrix Clustering ###
pdf('Run1_Zinbwave_Clustering.pdf')
data.frame(W.1, CellCyclePhase=colData(W.1)$Phase) %>%
  ggplot(aes(W1.1, W2.1, colour=CellCyclePhase)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()

pdf('Run2_Zinbwave_Clustering.pdf')
data.frame(W.2, CellCyclePhase=colData(W.2)$Phase) %>%
  ggplot(aes(W1.2, W2.2, colour=CellCyclePhase)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()

### Visualize Zinbwave tSNE Representation ###
set.seed(93024)

library(Rtsne)
r1_tsne_data <- Rtsne(W.1, pca = FALSE, perplexity=10, max_iter=5000)

pdf('Run1_Zinbwave_tSNE_Perplexity=10.pdf')
data.frame(Dim1=r1_tsne_data$Y[,1], Dim2=r1_tsne_data$Y[,2], 
           CellCyclePhase=ColData(W.1)$Phase) %>%
  ggplot(aes(Dim1, Dim2, colour=CellCyclePhase)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()

r1_tsne_data <- Rtsne(W.1, pca = FALSE, perplexity=30, max_iter=5000)

pdf('Run1_Zinbwave_tSNE_Perplexity=30.pdf')
data.frame(Dim1=r1_tsne_data$Y[,1], Dim2=r1_tsne_data$Y[,2], 
           CellCyclePhase=ColData(W.1)$Phase) %>%
  ggplot(aes(Dim1, Dim2, colour=CellCyclePhase)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()

r1_tsne_data <- Rtsne(W.1, pca = FALSE, perplexity=50, max_iter=5000)

pdf('Run1_Zinbwave_tSNE_Perplexity=50.pdf')
data.frame(Dim1=r1_tsne_data$Y[,1], Dim2=r1_tsne_data$Y[,2], 
           CellCyclePhase=ColData(W.1)$Phase) %>%
  ggplot(aes(Dim1, Dim2, colour=CellCyclePhase)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()

r2_tsne_data <- Rtsne(W.2, pca = FALSE, perplexity=10, max_iter=5000)

pdf('Run2_Zinbwave_tSNE_Perplexity=10.pdf')
data.frame(Dim1=r2_tsne_data$Y[,1], Dim2=r2_tsne_data$Y[,2], 
           CellCyclePhase=ColData(W.2)$Phase) %>%
  ggplot(aes(Dim1, Dim2, colour=CellCyclePhase)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()

r2_tsne_data <- Rtsne(W.2, pca = FALSE, perplexity=30, max_iter=5000)

pdf('Run2_Zinbwave_tSNE_Perplexity=30.pdf')
data.frame(Dim1=r2_tsne_data$Y[,1], Dim2=r2_tsne_data$Y[,2], 
           CellCyclePhase=ColData(W.2)$Phase) %>%
  ggplot(aes(Dim1, Dim2, colour=CellCyclePhase)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()

r2_tsne_data <- Rtsne(W.2, pca = FALSE, perplexity=50, max_iter=5000)

pdf('Run2_Zinbwave_tSNE_Perplexity=50.pdf')
data.frame(Dim1=r2_tsne_data$Y[,1], Dim2=r2_tsne_data$Y[,2], 
           CellCyclePhase=ColData(W.2)$Phase) %>%
  ggplot(aes(Dim1, Dim2, colour=CellCyclePhase)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()


### Find differential expression among cells in each run with DESeq2 ###
library(DESeq2)
library(apeglm)

### Run 1:
dds <- DESeqDataSetFromMatrix(countData = assay(r1_surf), colData = metadata(r1_surf), design = ~ seurat_clusters + Phase)
dds <- DESeq(dds, sfType="poscounts", useT = TRUE, minmu=1e-6)
res <- lfcShrink(dds, coef=c("Phase_G2M_vs_G1","Phase_S_vs_G1"), type = "apeglm")
head(res)

### Run 2:
dds <- DESeqDataSet(countData = assay(r2_surf), colData = metadata(r2_surf), design = ~ seurat_clusters + Phase)
dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
res <- lfcShrink(dds, coef=c("Phase_G2M_vs_G1", "Phase_S_vs_G1"), type = "apeglm")
head(res)



### Find integration anchors in the datasets and integrate the two runs together ###
df_anchors <- FindIntegrationAnchors(object.list=c(r1,r2), dims=1:30)
df <- IntegrateData(anchorset= df_anchors, dims=1:30)

### Apply the same Zinbwave algorithm to the integrated dataset ###
SingleCellExperiment(df[["RNA"]]@data)->df.sce
filter.int <- rowSums(assay(df.sce)>5)>5
table(filter.int)
# filter.int
df.sce <- df.sce[filter.int,]
rowVars(as.matrix(assay(df.sce))) -> vars.int
names(vars.int) <- rownames(df.sce)
vars.int <- sort(vars.int, decreasing = TRUE)
head(vars.int)
df.sce <- df.sce[names(vars.int)[1:2000],]
assayNames(df.sce)[1] <- "counts"
as.matrix(assay(df.sce))->assay(df.sce)
ifelse(assay(df.sce)==".",0,assay(df.sce))->assay(df.sce)
df_zinb <- zinbwave(df.sce, K = 20, epsilon = 1000000, normalizedValues = T, residuals = T, verbose = T)
W.int <- reducedDim(df_zinb)
df_zinb
write.table(W.int,"BothRuns_ZinbWave_W.xlsx",quote=F,sep="\t",col.names=T,row.names=T)

### Visualize the W Matrix Clustering ###
pdf('Combined_Zinbwave_Clustering.pdf')
data.frame(W.int, CellCyclePhase=colData(W.int)$Phase) %>%
  ggplot(aes(W1.int, W2.int, colour=CellCyclePhase)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()

### Visualize ZinbWave tSNE Representation ###
df_tsne_data <- Rtsne(W.int, pca = FALSE, perplexity=10, max_iter=5000)

pdf('Combined_Zinbwave_tSNE_Perplexity=10.pdf')
data.frame(Dim1=df_tsne_data$Y[,1], Dim2=df_tsne_data$Y[,2], 
           CellCyclePhase=ColData(W.int)$Phase) %>%
  ggplot(aes(Dim1, Dim2, colour=CellCyclePhase)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()

df_tsne_data <- Rtsne(W.int, pca = FALSE, perplexity=30, max_iter=5000)

pdf('Combined_Zinbwave_tSNE_Perplexity=30.pdf')
data.frame(Dim1=df_tsne_data$Y[,1], Dim2=df_tsne_data$Y[,2], 
           CellCyclePhase=ColData(W.int)$Phase) %>%
  ggplot(aes(Dim1, Dim2, colour=CellCyclePhase)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()

df_tsne_data <- Rtsne(W.int, pca = FALSE, perplexity=50, max_iter=5000)

pdf('Combined_Zinbwave_tSNE_Perplexity=50.pdf')
data.frame(Dim1=df_tsne_data$Y[,1], Dim2=df_tsne_data$Y[,2], 
           CellCyclePhase=ColData(W.int)$Phase) %>%
  ggplot(aes(Dim1, Dim2, colour=CellCyclePhase)) + geom_point() + 
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()


dds <- DESeqDataSet(countData = assay(df_surf), colData = metadata(df_surf), design = ~ seurat_clusters + Phase)
dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
res <- lfcShrink(dds, coef=c("Phase_G2M_vs_G1", "Phase_S_vs_G1"), type = "apeglm")
head(res)

### Save SingleCellExperiment object with Zinbsurf reduction as a Seurat Object and Export as RDS ###
as.seurat(r1_surf, counts = "counts", data ="counts")-> r1_seu
r1_seu <- FindNeighbors(r1_seu, reduction = "zinbwave", dims = 1:2)
seu <- FindClusters(object = r1_seu)
saveRDS(file = "Run1_ZinbSurf_Top1000_Epsilon=1000000_CellCyclePhase_DE.rds")

as.seurat(r2_surf, counts = "counts", data = "counts")->r2_seu
r2_seu <- FindNeighbors(r2_seu, reduction ="Zinbwave", dims = 1:2)
r2_seu <- FindClusters(object = r2_seu)
saveRDS(r2_seu, file = "Run2_ZinbSurf_Top1000_Epsilon=1000000_CellCyclePhase_DE.rds")

as.seurat(df_surf, counts = "counts", data ="counts")-> df_seu
df_seu <- FindNeighbors(df_seu, reduction = "Zinbwave", dims = 1:2)
df_seu <- FindClusters(object = df_seu)
saveRDS(df_seu, file = "Combined_ZinbSurf_Top1000_Epsilon=1000000_CellCyclePhase_DE.rds")


### Setup data for iDEA Analysis for Gene Set Enrichment of Single Cell Data ###
library(iDEA)
data("mouseGeneSets")
as.data.frame(res[,c(2,3)])-> Sum


