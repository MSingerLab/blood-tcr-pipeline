#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(viridis)
library(gridExtra)
library(RColorBrewer)
#1: 10X dir, 2: freq file, 3: sample, 4: output folder loc
#Example: Rscript pipeline.R /Users/jacobluber/Desktop/ucsf/working_soLNGEX/outs/filtered_gene_bc_matrices/GRCh38 /Users/jacobluber/Dropbox\ \(HMS\)/Jacob\ Work/Analysis/layn/chi-square/working_so_plotting.csv working_so /Users/jacobluber/Dropbox\ \(HMS\)/Jacob\ Work/working/working_so/
args = commandArgs(trailingOnly=TRUE)
tenxdir <- args[1]
freq <- args[2]
outloc <- args[3]
sample <- args[4]

#optionally read in existing Seurat object
working_so <- readRDS(tenxdir)

#optionally use only certain clusters and regenerate UMAP
#working_so <- SubsetData(working_so, ident.use = c(0,2,14,9,7,5,15,6,10,4,13), do.clean = TRUE, do.scale = TRUE)
#all.genes <- rownames(x = working_so)
#working_so <- ScaleData(object = working_so, features = all.genes)
#working_so <- FindVariableFeatures(object = working_so, mean.function = ExpMean, dispersion.function = LogVMR,  do.plot = FALSE)
#hv.genes <- head(x = VariableFeatures(object = working_so), 1000)
#working_so <- RunPCA(object = working_so, pc.genes = hv.genes, do.print = FALSE, pcs.print = 1:5,genes.print = 5, pcs.compute = 50)
#working_so <- FindNeighbors(object = working_so, dims = 1:30)
#working_so <- FindClusters(object = working_so, resolution = 1.2)
#working_so <- RunUMAP(object = working_so, reduction.use = "pca", dims = 1:15, n_neighbors = 15, min_dist = 0.3)



tcr <- read.csv(freq)
tcr.idents <- tcr$Group
tcr$Barcode <- gsub("(.*)\\-(1)","\\1",tcr$Barcode)
names(tcr.idents)=tcr$Barcode
working_so <- AddMetaData(working_so,tcr.idents,col.name='matching')
working_so <- AddMetaData(working_so,data.frame(working_so@reductions$umap@cell.embeddings)$UMAP_1,col.name='u1')
working_so <- AddMetaData(working_so,data.frame(working_so@reductions$umap@cell.embeddings)$UMAP_2,col.name='u2')
new.idents <- tcr$Tcr
names(new.idents)=tcr$Barcode
working_so <- AddMetaData(working_so,new.idents,col.name="TCR")
new.idents <- tcr$Frequency
names(new.idents)=tcr$Barcode
working_so <- AddMetaData(working_so,new.idents,col.name="Frequency")
#vals <- working_so@meta.data$Clone
#replacevals <- replace(vals,vals==0,NA)
#working_so@meta.data$Clone <- replacevals


saveRDS(working_so, file = paste(outloc,sample,".rds",sep=''))



