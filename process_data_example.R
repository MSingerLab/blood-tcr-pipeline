#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(viridis)
library(gridExtra)
library(RColorBrewer)
library(tibble)
#1: 10X dir, 2: freq file, 3: sample, 4: output folder loc
#Example: Rscript pipeline.R /Users/jacobluber/Desktop/ucsf/working_soLNGEX/outs/filtered_gene_bc_matrices/GRCh38 /Users/jacobluber/Dropbox\ \(HMS\)/Jacob\ Work/Analysis/layn/chi-square/working_so_plotting.csv working_so /Users/jacobluber/Dropbox\ \(HMS\)/Jacob\ Work/working/working_so/
#parse command line data
#args = commandArgs(trailingOnly=TRUE)
#tenxdir <- args[1]
#freq <- args[2]
#sample <- args[3]
setwd("/Users/jacobluber/Desktop/cwd/k468")

data <- Read10X(data.dir = "/Volumes/Jacob's Backup /k468")
#read data from Kraken and create Seurat object
#data <- Read10X(data.dir = tenxdir)

#data <- Read10X(data.dir = "/Volumes/Jacob's Backup /k468")
data2 <- CreateSeuratObject(data, min.cells = 3,min.features=400, project = "test")
housekeeping_genes_pre <- read.table("/Users/jacobluber/Desktop/cwd/HK_Satija.txt",stringsAsFactors = FALSE)
housekeeping_genes <- housekeeping_genes_pre$V1
housekeeping_genes <- intersect(rownames(data2), housekeeping_genes)
subset_counts_matrix <- data.frame(data2@assays$RNA@counts)[housekeeping_genes,]
subset_counts_matrix[subset_counts_matrix > 0] =1 
summary(colSums(subset_counts_matrix))
hist(colSums(subset_counts_matrix),main="> .5 selected as cutoff for HK genes")
pdf("HK_genes_hist.pdf",width=7,height=7)
hist(colSums(subset_counts_matrix),main="> .5 selected as cutoff for HK genes")
dev.off()
counts <- colSums(subset_counts_matrix)
counts2 <- data.frame(counts)
counts3 <- counts2 %>% tibble::rownames_to_column() %>% filter(counts > 48)


all.genes <- rownames(x = data2)
mito_pre <- read.table("/Users/jacobluber/Desktop/mito_genes.txt",stringsAsFactors = FALSE)
mito_genes <- mito_pre$V1
mito_genes <- intersect(rownames(data2), mito_genes)
subset_counts_matrix <- data.frame(data2@assays$RNA@counts)[mito_genes,]
subset_counts_matrix[subset_counts_matrix > 0] =1 
summary(colSums(subset_counts_matrix))
hist(colSums(subset_counts_matrix),main="< 200/1157 genes (Broad list) selected as cutoff for mito genes")
pdf("mito_genes_hist.pdf",width=7,height=7)
hist(colSums(subset_counts_matrix),main="< 200/1157 genes (Broad list) selected as cutoff for mito genes")
dev.off()

counts4 <- colSums(subset_counts_matrix)
counts5 <- data.frame(counts4)
counts6 <- counts5 %>% tibble::rownames_to_column() %>% filter(counts4 < 200)

cells <- intersect(counts6$rowname,counts3$rowname)

#Run all processing through UMAP step
#all.genes <- rownames(x = working_so)
#housekeeping_genes_pre <- read.table("/Users/jacobluber/Desktop/cwd/HK_Satija.txt")
#housekeeping_genes <- housekeeping_genes_pre$V1
#non.mito.genes <- grep("^MT-", x = all.genes, value = TRUE, invert = TRUE)
#mito.genes <- grep("^MT-", x = all.genes, value = TRUE, invert = FALSE)

#housekkeeping and mito gene plots
# hist_arr_hk = c()
# for (i in 1:length(colnames(data2@assays$RNA@counts))) {
#     hk_genes <- 0
#     for (j in 1:length(housekeeping_genes)) {
#       if (data2@assays$RNA@counts[,i][housekeeping_genes[j]] == 1){
#         old <- hk_genes
#         hk_genes <- old+1
#       }
#     }
#     hist_arr_hk<- c(hist_arr_hk,hk_genes)
# }

# hist_arr_hk2 = c()
# for (i in 1:100) {
#   print(i)
#   hk_genes <- 0
#   for (j in 1:length(housekeeping_genes)) {
#     if (data2@assays$RNA@counts[,i][housekeeping_genes[j]] >= 1){
#       hk_genes <- hk_genes+1
#     }
#   }
#   hist_arr_hk2<- c(hist_arr_hk2,hk_genes/106)
# }
# 
# m_genes_arr <- c()
# for (i in 1:300) {
#   print(i)
#   m_genes <- 0
#   for (j in 1:length(mito.genes)) {
#     if (data2@assays$RNA@counts[,i][mito.genes[j]] == 1){
#       m_genes <- m_genes+1
#     }
#   }
#   m_genes_arr<- c(m_genes_arr,m_genes/13)
# }

working_so <- NormalizeData(data2)
working_so <- SubsetData(working_so, cells = cells)
working_so <- ScaleData(object = working_so, features = rownames(working_so))
working_so <- FindVariableFeatures(object = working_so, mean.function = ExpMean, dispersion.function = LogVMR,
                                   do.plot = FALSE)
hv.genes <- head(x = VariableFeatures(object = working_so), 1000)
working_so <- RunPCA(object = working_so, pc.genes = hv.genes, do.print = FALSE, pcs.print = 1:5,
                     genes.print = 5, pcs.compute = 50)
working_so <- FindNeighbors(object = working_so, dims = 1:30)
working_so <- FindClusters(object = working_so, resolution = 1.2)
working_so <- RunUMAP(object = working_so, reduction.use = "pca", dims = 1:15, n_neighbors = 15, min_dist = 0.3)

#optionally read in existing Seurat object
#working_so <- readRDS("/Users/jacobluber/browser/k409.ln-blood.rds")

#get the tcr data to determine matching cluster status, and add various other things to the metadata 
tcr <- read.csv("/Users/jacobluber/browser/k468.blood-ln.csv")
#tcr <- read.csv(freq)
tcr.idents <- tcr$matching
#CHECK INPUT FILE
tcr$pcluster <- as.numeric(tcr$matching)
tcr$pcluster[tcr$pcluster==1] <- 1
tcr$pcluster[tcr$pcluster==2] <- 1
tcr$pcluster[tcr$pcluster==3] <- 1
tcr$pcluster[tcr$pcluster==4] <- 2
tcr$pcluster[tcr$pcluster==5] <- 2
tcr$ambig <- as.numeric(tcr$matching)
tcr$ambig[tcr$ambig==1] <- 1
tcr$ambig[tcr$ambig==2] <- 1
tcr$ambig[tcr$ambig==3] <- 2
tcr$ambig[tcr$ambig==4] <- 2
tcr$ambig[tcr$ambig==5] <- 2
tcr$consis <- as.numeric(tcr$matching)
tcr$consis[tcr$consis==1] <- 2
tcr$consis[tcr$consis==2] <- 2
tcr$consis[tcr$consis==3] <- 1
tcr$consis[tcr$consis==4] <- 2
tcr$consis[tcr$consis==5] <- 2
tcr$Barcode <- tcr$`X`
tcr$Barcode <- gsub("(.*)\\-(1)","\\1",tcr$Barcode)
names(tcr.idents)=tcr$Barcode
tcr.pclusters <- tcr$pcluster
tcr.ambig <- tcr$ambig
names(tcr.pclusters)=tcr$Barcode
names(tcr.ambig)=tcr$Barcode
tcr.consis <- tcr$consis
names(tcr.consis)=tcr$Barcode
working_so <- AddMetaData(working_so,tcr.pclusters,col.name='pcluster')
working_so <- AddMetaData(working_so,tcr.consis,col.name='consistent')
working_so <- AddMetaData(working_so,tcr.ambig,col.name='ambig')
working_so <- AddMetaData(working_so,tcr.idents,col.name='matching')
working_so <- AddMetaData(working_so,data.frame(working_so@reductions$umap@cell.embeddings)$UMAP_1,col.name='u1')
working_so <- AddMetaData(working_so,data.frame(working_so@reductions$umap@cell.embeddings)$UMAP_2,col.name='u2')


#determine which clusters to keep 
clusters <- levels(working_so@meta.data$seurat_clusters)

#iterative function to see if a given cluster should be kept, results aggregate like a fold in Haskell 
see_if_keep_cluster <- function(so, cluster, idents){
  copy_so <- so
  cso <- SubsetData(copy_so, ident.use = cluster, do.clean = TRUE, do.scale = TRUE)
  new_id <- NULL
  cd3e <- sum(GetAssayData(object = cso, slot = "data")["CD3E",]>0)/nrow(cso@meta.data)
  cd3d <- sum(GetAssayData(object = cso, slot = "data")["CD3D",]>0)/nrow(cso@meta.data)
  cd3g <- sum(GetAssayData(object = cso, slot = "data")["CD3G",]>0)/nrow(cso@meta.data)
  cd8b <- sum(GetAssayData(object = cso, slot = "data")["CD8B",]>0)/nrow(cso@meta.data)
  cd8a <- sum(GetAssayData(object = cso, slot = "data")["CD8A",]>0)/nrow(cso@meta.data)
  foxp3 <- sum(GetAssayData(object = cso, slot = "data")["FOXP3",]>0)/nrow(cso@meta.data)
  cd4 <- sum(GetAssayData(object = cso, slot = "data")["CD4",]>0)/nrow(cso@meta.data)
  if (cd3e > .5 || cd3d > .5 || cd3g > .5) {
    if (cd8b > .5 && cd8a > .5 && foxp3 < .05 && cd4 < .05) {
      new_id <- cluster 
    }
  }
  idents <- c(idents,new_id)
  return(idents)
}

idents <- c()
for (i in clusters){
  idents <- see_if_keep_cluster(working_so, as.numeric(i),idents)
}

#CHANGE
sample <- "k468"

working_so1 <- SubsetData(working_so, ident.use = idents, do.clean = TRUE, do.scale = TRUE)
working_so1 <- ScaleData(object = working_so1, features = rownames(working_so1))
working_so1 <- FindVariableFeatures(object = working_so1, mean.function = ExpMean, dispersion.function = LogVMR,
                                    do.plot = FALSE)
hv.genes <- head(x = VariableFeatures(object = working_so1), 1000)
working_so1 <- RunPCA(object = working_so1, pc.genes = hv.genes, do.print = FALSE, pcs.print = 1:5,
                      genes.print = 5, pcs.compute = 50)
working_so1 <- FindNeighbors(object = working_so1, dims = 1:30)
working_so1 <- FindClusters(object = working_so1, resolution = 1.2)
working_so1 <- RunUMAP(object = working_so1, reduction.use = "pca", dims = 1:15, n_neighbors = 15, min_dist = 0.3)

#Generate UMAPs to check cluster selection 
t1 <- paste0("Clusters Kept=",paste(idents, collapse=', ' ))
p1 <- DimPlot(object = working_so, reduction.use = "umap", no.legend = FALSE, do.return = TRUE, label = TRUE, vector.friendly = TRUE, pt.size = .4) + ggtitle(paste0(length(colnames(data2))," -> ",length(colnames(working_so))," cells (HK/mito filtering)")) + theme(plot.title = element_text(hjust = 0.5)) + coord_fixed() 
p2 <- DimPlot(object = working_so1, reduction.use = "umap", no.legend = FALSE, do.return = TRUE, label = TRUE, vector.friendly = TRUE, pt.size = .4) + ggtitle(paste0(length(colnames(working_so))," -> ",length(colnames(working_so1))," cells (CD8 filtering)")) + theme(plot.title = element_text(hjust = 0.5)) + coord_fixed() + labs(subtitle = t1)
p3 <- FeaturePlot(working_so1,features=c("FOXP3","CD4","CD3G","CD8B","CD8A","CD3D","CD3E"), reduction = "umap",cols=c("gray","red"),pt.size = .4, coord.fixed = TRUE,combine=FALSE)
title <- paste0(sample,"_UMAP_original.pdf")
pdf(title,width=7,height=7)
p1
dev.off()
pdf(paste0(sample,"_UMAP_cd8s.pdf"),width=7,height=7)
p2
dev.off()
pdf(paste0(sample,"_UMAP_genes.pdf"),width=7,height=7)
p3
dev.off()
pdf(paste0(sample,"_UMAP_matching_plots.pdf"),width=7,height=7)
FeaturePlot(working_so1,features=c("pcluster"), reduction = "umap",pt.size = .4, coord.fixed = TRUE) + NoLegend() + labs(title = "Definition Used (consistent U ambiguous), matching in grey")
FeaturePlot(working_so1,features=c("ambig"), reduction = "umap",pt.size = .4, coord.fixed = TRUE) + NoLegend() + labs(title = "Ambiguous Definition, matching in grey")
FeaturePlot(working_so1,features=c("consistent"), reduction = "umap",pt.size = .4, coord.fixed = TRUE) + NoLegend() + labs(title = "Consistent Definition, matching in grey")
dev.off()

#rerun UMAP pipeline on data with non CD8 clusters removed 
#working_so1 <- SubsetData(working_so, ident.use = idents, do.clean = TRUE, do.scale = TRUE)


markers <- FindMarkers(working_so1, ident.1 = 1, ident.2 = NULL, only.pos = TRUE,test.use = "MAST")
write.table(markers,"MAST_matching_markers.txt")

w#write output files for COMET 
write.table(GetAssayData(object = working_so1, slot = "counts"),file=paste0("counts_matrix",".txt"),sep="\t",quote=FALSE,row.names=TRUE)
write.table(GetAssayData(object = working_so1, slot = "scale.data"),file=paste0("normalized_matrix",".txt"),sep="\t",quote=FALSE,row.names=TRUE)
write.table(Embeddings(working_so1[["umap"]]),file=paste0("viz",".txt"),sep="\t",quote=FALSE,row.names=TRUE)
cluster_out <- working_so1@meta.data %>% tibble::rownames_to_column('barcode') %>% select(barcode,pcluster)
write.table(cluster_out,file=paste0("clusters",".txt"),sep="\t",quote=FALSE,row.names=FALSE)
write.csv(working_so1@meta.data,file=paste0("metadata",".csv"))
#write out metadata
saveRDS(working_so, file = paste0(sample,".rds"))

#p1 <- DimPlot(object = working_so, reduction.use = "umap", no.legend = FALSE, do.return = TRUE, label = TRUE, vector.friendly = TRUE, pt.size = 1) + ggtitle("UMAP Clusters Removed") + theme(plot.title = element_text(hjust = 0.5)) + coord_fixed()
p4 <- FeaturePlot(working_so,features=c("FOXP3","CD4","CD3G","CD8B","CD8A","CD3D","CD3E"), reduction = "umap",cols=c("gray","red"),pt.size = .4, coord.fixed = TRUE,combine=FALSE)
title1 <- paste0(sample,"_UMAP_genes_noCD8_filtering.pdf")
pdf(title1,width=20,height=6)
p4
dev.off()
