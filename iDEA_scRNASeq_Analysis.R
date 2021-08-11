### Prepare the SummaryData and AnnotatedData files for iDEA ###

library(iDEA)
library(dplyr)
library(tidyverse)
library(biomaRt)

read.table("Run1/Zinbwave/DiffExp/Run1_ZinbWave_Top2000_Epsilon=1000000_Cluster0vs1_DE.txt", header = T, sep = "\t")->Sum
data.frame(mgi_symbol = Sum$mgi_symbol, beta = Sum$log2FoldChange, beta_var = Sum$lfcSE)->Sum
head(Sum)

data(mouseGeneSets)
mouseGeneSets[1:5,1:5]
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
genes <- getBM(filter = "mgi_symbol", attributes = c("mgi_symbol","mgi_id"), mart = mart, values = Sum$mgi_symbol, verbose = T)
data <- merge(Sum, genes, by.x = "mgi_symbol", by.y = "mgi_symbol")
data.frame(beta = data$beta, beta_var = data$beta_var, mgi_id = data$mgi_id)->data
rownames(data)<-data$mgi_id
data[,1:2] -> data
head(data)
write.table(data, "iDEA_SummaryData.txt", row.names = T, col.names = T, sep = "\t", quote = F)

### Create iDEA Object and Run Analysis ###
data->summary
mouseGeneSets->annotation
max_var_beta = 100
min_precent_annot = 0.0025
num_core = 1
if(!is.data.frame(summary)){summary <- as.data.frame(summary)}
print('Checking data format')

colnames(summary) <- c("beta", "beta_var")
keep_index <- summary$beta_var<max_var_beta & !is.na(summary$beta_var)
summary <- as.data.frame(summary[keep_index,])
if(!is.data.frame(annotation)){annotation <- as.data.frame(annotation)}
print('Data format correct')

if(nrow(summary) != nrow(annotation)){
  summary$GeneID <- rownames(summary)
  annotation$GeneID <- rownames(annotation)
  combined_data <- merge(summary, annotation, by="GeneID")
  summary <- combined_data[, 2:3]
  rownames(summary) <- combined_data$GeneID
  annotation <- combined_data[, 4:ncol(combined_data)]
  rownames(annotation) <- combined_data$GeneID
  rm(combined_data)}
print('Data merged')

num_gene <- nrow(annotation)
precent_annot <- apply(annotation, 2, sum)/num_gene
print('Data filtered')

annotation <- annotation[, which(precent_annot>min_precent_annot)]
annotation = as.data.frame(annotation)
annot_list <- lapply(1:ncol(annotation), FUN = function(x) { which(annotation[,x]==1) })
print('Data converted to list')

names(annot_list) <- colnames(annotation)
idea <- new(
  Class = "iDEA",
  gene_id = rownames(annotation),
  annot_id = colnames(annotation),
  summary = summary,
  annotation = annot_list,
  num_gene = num_gene)
print('iDEA Object Created')
idea

#idea <- CreateiDEAObject(Sum, mouseGeneSets, max_var_beta = 100, min_precent_annot = 0.0025, num_core = 1)
head(idea@summary)
head(idea@annotation[[1]])

idea@num_core <- as.numeric(1)
idea <- iDEA.fit(idea, fit_noGS=FALSE, init_beta=NULL, init_tau=c(-2,0.5), min_degene=5, em_iter=15, mcmc_iter=1000, fit.tol=1e-5, modelVariant = F, verbose=TRUE)

### Now do the iDEA fit algorithm to fit the iDEA model onto the input data ###
# Set parameters: fit_noGS=FALSE, init_beta=NULL, modelVariant = F #
init_tau=c(-2,0.5)
min_degene=5
em_iter=15
mcmc_iter=1000
fit.tol=1e-5
idea@num_core <- as.numeric(1)

num_core <- idea@num_core
num_gene <- idea@num_gene
num_annot <- length(idea@annotation)

cat(paste("## ===== iDEA INPUT SUMMARY ==== ##\n"))
cat(paste("## number of annotations: ", num_annot,"\n"))
cat(paste("## number of genes: ", num_gene,"\n"))
init_beta <- idea@summary[,1]
Annot <- as.matrix(data.frame(rep(1, num_gene)))
res_idea <- lapply(1:num_annot, FUN = function(x) {Annot <- rep(0, idea@num_gene)
Annot[idea@annotation[[x]]] <- 1
Annot <- Annot - mean(Annot)
Annot <- as.matrix(data.frame(rep(1, num_gene), Annot))
t1 <- system.time(model1 <- try(res <- EMMCMCStepSummary(idea@summary[,1], idea@summary[,2], as.matrix(Annot), init_beta, init_tau, em_iter, mcmc_iter, min_degene)))
if(class(model1) != "try-error"){
  rownames(res$pip) <- idea@gene_id
  colnames(res$pip) <- "PIP"
  rownames(res$beta) <- idea@gene_id
  colnames(res$beta) <- "beta"
  rownames(res$annot_coef) <- c("tau_1", "tau_2")
  colnames(res$annot_coef) <- "annot_coef"
  rownames(res$annot_var) <- c("tau_1","tau_2")
  colnames(res$annot_var) <- "annot_var"
  res$converged <- TRUE
  res$cttime <- t1[3]}
print('The loop is completed')
}
)

res_idea <- lapply(1:num_annot, FUN = function(x) {
Annot <- rep(0, idea@num_gene)
Annot[idea@annotation[[x]]] <- 1
Annot <- Annot - mean(Annot)
Annot <- as.matrix(data.frame(rep(1, num_gene), Annot) )
t1 <- system.time(model1 <- try( res <- EMMCMCStepSummary(idea@summary[,1], idea@summary[,2], as.matrix(Annot), init_beta, init_tau, em_iter, mcmc_iter, min_degene) ))
if(class(model1) != "try-error"){
  rownames(res$pip) <- idea@gene_id
  colnames(res$pip) <- "PIP"
  rownames(res$beta) <- idea@gene_id
  colnames(res$beta) <- "beta"
  rownames(res$annot_coef) <- c("tau_1", "tau_2")
  colnames(res$annot_coef) <- "annot_coef"
  rownames(res$annot_var) <- c("tau_1", "tau_2")
  colnames(res$annot_var) <- "annot_var"
  res$converged   <- TRUE
  res$ctime   <- t1[3]
}
}
)

#t1 <- system.time(model1 <- try( res <- EMMCMCStepSummaryVariant(idea@summary[,1], idea@summary[,2], as.matrix(Annot), init_beta, init_tau[1], em_iter, mcmc_iter, min_degene) ))
#                 if(class(model1) != "try-error"){
#                     rownames(res$pip) <- idea@gene_id
#                     res$converged   <- TRUE
#                     res$ctime   <- t1[3]
#                 }else{res <- NULL}
#                 idea@noGS <- res
#
names(res_idea) <- idea@annot_id
idea@de <- res_idea
rm(res_idea)

idea

head(idea@de)
idea <- iDEA.louis(idea)
head(idea@gsea)

write.table(idea@de, "iDEA_GSEA_Analysis.txt", quote = F, row.names = T, col.names = T, sep = "\t")
  