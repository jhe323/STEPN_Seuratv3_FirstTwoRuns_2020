library(patchwork)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(dplyr)
library(cowplot)

workspace<-"/home/CAM/jhe/scRNASeq/Analysis/Combined/"

setwd(workspace)
filt <- readRDS("./Combined_FUS1_scRNASeq_Seurat.rds")
#pdf('MAP3K7_Vplot.pdf')
#VlnPlot(filt,features=c("Map3k7"),pt.size=0)
#pdf('SMAD2_Vplot.pdf')
#VlnPlot(filt,features=c("Smad2"),pt.size=0)
#pdf('SMAD3_Vplot.pdf')
#VlnPlot(filt,features=c("Smad3"),pt.size=0)
#pdf('SMAD7_Vplot.pdf')
#VlnPlot(filt,features=c("Smad7"),pt.size=0)
pdf('SMAD4_Vplot.pdf')
VlnPlot(filt,features=c("Smad4"),pt.size=0)

#pdf('SMAD2_Fplot.pdf')
#FeaturePlot(filt,features=c("Smad2"),pt.size=0.25)
#pdf('SMAD3_Fplot.pdf')
#FeaturePlot(filt,features=c("Smad3"),pt.size=0.25)
#pdf('SMAD7_Fplot.pdf')
#FeaturePlot(filt,features=c("Smad7"),pt.size=0.25)
#pdf('MAP3K7_Fplot.pdf')
#FeaturePlot(filt,features=c("Map3k7"),pt.size=0.25)
pdf('SMAD4_Fplot.pdf')
FeaturePlot(filt,features=c("Smad4"),pt.size=0.25)
