#!/usr/bin/env Rscript

# README ------------------------------------------------------------------
# Assign some cells in Aza-treated co-culture as SNU719 or NCC24,
# using the monoculture Aza-treated SNU719 and monoculture Aza-treated NCC24 cells as a reference.

# Set-up packages, working directory, datasets -----------------------------------------------------------
## Load Packages
library(Seurat)
library(sctransform)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(cowplot)
library(Matrix)
library(scales)
library(bitops)
library(RCurl)
library(Matrix.utils)
library(SingleCellExperiment)
library(beepr)

options(future.globals.maxSize = 4000 * 1024^2)

## Set working directory
setwd("/scrnaseq/results/cell line assignment/Aza")

## Load datasets
# Load Aza-treated SNU719 in monoculture and Aza-treated NCC24 in monoculture for building reference dataset
load("/scrnaseq/results/dataset preparation/mono/april2017/SNU719-Aza-april2017/SCT-SNU719-Aza-april2017.RData")
load("/scrnaseq/results/dataset preparation/mono/sept2017/SNU719-Aza-sept2017/SCT-SNU719-Aza-sept2017.RData")
load("/scrnaseq/results/dataset preparation/mono/sept2017/NCC24-DMSO-sept2017/SCT-NCC24-DMSO-sept2017.RData")
# Load DMSO-treated co-culture
load("/scrnaseq/results/dataset preparation/co/co-Aza/SCT-co_Aza.RData")
# Forming reference dataset -----------------------------------------------------------
# Merge SNU719 and NCC24  to form a reference (object: monoculture.Aza.merged)
monoculture.Aza.merged <- merge(x = NCC24.DMSO.sept2017.seuratobject,
                         y = c(SNU719.Aza.april2017.seuratobject,SNU719.Aza.sept2017.seuratobject),
                         add.cell.id = c("NCC24.DMSO.sept2017","SNU719.Aza.april2017","SNU719.Aza.sept2017"),
                         project = "merge_ref_Aza")
# Run SCTransform
monoculture.Aza.merged <-SCTransform(monoculture.Aza.merged,
                                      return.only.var.genes = F,
                                      vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
# Run PCA and plot Elbow Plot
monoculture.Aza.merged <- RunPCA(monoculture.Aza.merged, verbose = FALSE)
png("Elbow Plot of merged Aza-treated SNU719 and NCC24 in monoculture.png")
ElbowPlot(monoculture.Aza.merged) + ggtitle("Merged Aza-treated SNU719 and NCC24 in monoculture")
dev.off()

png("PCA Plot of merged Aza-treated SNU719 and NCC24 in monoculture.png")
PCAPlot(monoculture.Aza.merged, group.by = "sample") +
  ggtitle("Merged Aza-treated SNU719 and NCC24 in monoculture")
dev.off()

monoculture.Aza.merged <- RunUMAP(monoculture.Aza.merged, reduction = "pca",dims = 1:15)
monoculture.Aza.merged <- FindNeighbors(monoculture.Aza.merged, dims = 1:15)
monoculture.Aza.merged <- FindClusters(monoculture.Aza.merged, resolution = c(0.1))

# Observe clustering with UMAP at resolution of 0.1
Idents(object = monoculture.Aza.merged) <- "SCT_snn_res.0.1"
png("UMAP of merged Aza-treated SNU719 and NCC24 in monoculture.png")
DimPlot(monoculture.Aza.merged, reduction = "umap", group.by = "CellLine",
        cols = c('SNU719' = 'tomato2','NCC24' = 'skyblue2')) +
  ggtitle("Merged Aza-treated SNU719 and NCC24 in monoculture")
dev.off()

# Save reference dataset
save(monoculture.Aza.merged, file="monoculture_Aza_merged.RData")

# Cell line assignment in co-culture -----------------------------------------------------------
reference.dataset <- monoculture.Aza.merged
query.dataset <- coculture.Aza

DefaultAssay(query.dataset) <- "SCT"

# Set reduction to CCA for FindTransferAnchors
query.anchors <- FindTransferAnchors(reference = reference.dataset, query = query.dataset, dims = 1:10,
                                     reduction = 'cca')
predictions <- TransferData(anchorset = query.anchors, refdata = reference.dataset$CellLine, dims = 1:10,
                            weight.reduction = 'cca')
query.dataset <- AddMetaData(query.dataset, metadata = predictions)

# Save predicted.id results to a csv file
write.csv(query.dataset@meta.data[["predicted.id"]], file = "Predicted Cell Line Assignment Aza.csv")

### Label Top 1500 cells with highest prediction score for SNU719 and NCC24
# Transfer predictions and query dataset to new objects
predictions.sc1 <- predictions
query.dataset.sc1 <- query.dataset
# Create a column to indicate which cells to label
predictions.sc1$predicted.id.sc1 <- "Unassigned"
# Sort by ordered SNU719 prediction scores, and label top 1500 cells as SNU719
predictions.sc1 <-predictions.sc1[order(-predictions.sc1$prediction.score.SNU719), ]
predictions.sc1$predicted.id.sc1[1:1500] <- "SNU719"
# Sort by ordered NCC24 prediction scores, and abel top 1500 cells as NCC24
predictions.sc1 <-predictions.sc1[order(-predictions.sc1$prediction.score.NCC24), ]
predictions.sc1$predicted.id.sc1[1:1500] <- "NCC24"
## Add metadata back to query.dataset.sc1
query.dataset.sc1 <- AddMetaData(query.dataset.sc1, metadata = predictions.sc1)
## Visualize with UMAP
query.dataset.sc1 <- RunUMAP(query.dataset.sc1, reduction = "pca", dims = 1:15)
query.dataset.sc1 <- FindNeighbors(query.dataset.sc1, dims = 1:15)
query.dataset.sc1 <- FindClusters(query.dataset.sc1, resolution = c(0.5), dims = 1:15)
Idents(object = query.dataset.sc1) <- "SCT_snn_res.0.5"

png("UMAP of top 1500 cells likely to be SNU719 and NCC24 for Aza co-culture.png")
DimPlot(query.dataset.sc1, reduction = "umap", group.by = "predicted.id.sc1", label = "FALSE",
        cols = c('SNU719' = 'tomato2','NCC24' = 'skyblue2', 'Unassigned' = 'gray70')) +
  ggtitle("Top 1500 cells likely to be SNU719 and NCC24 for Aza co-culture")
dev.off()

png("UMAP of top 1500 cells likely to be SNU719 and NCC24 for Aza co-culture on same axis as monoculture.png")
DimPlot(query.dataset.sc1, reduction = "umap", group.by = "predicted.id.sc1", label = "FALSE",
        cols = c('SNU719' = 'tomato2','NCC24' = 'skyblue2', 'Unassigned' = 'gray70')) +
  xlim(-10,5) +
  ylim(-4,4) +
  ggtitle("Top 1500 cells likely to be SNU719 and NCC24 for Aza co-culture")
dev.off()

# Labeling cells in co-culture as SNU719 and NCC24------------------------------------------------------------------
query.dataset <- RunUMAP(query.dataset, reduction = "pca", dims = 1:15)
query.dataset <- FindNeighbors(query.dataset, dims = 1:15)
query.dataset <- FindClusters(query.dataset, resolution = c(0.5), dims = 1:15)
# Observe clustering at resolution of 0.4
Idents(object = query.dataset) <- "SCT_snn_res.0.5"

png("UMAP of cells likely to be SNU719 and NCC24 for Aza co-culture.png")
DimPlot(query.dataset, reduction = "umap", group.by = "predicted.id", label = "TRUE",
        cols = c('SNU719' = 'tomato2','NCC24' = 'skyblue2')) +
  ggtitle("If all cells are SNU719 or NCC24 for Aza co-culture")
dev.off()

png("UMAP of unlabeled Aza co-culture.png")
DimPlot(query.dataset, reduction = "umap", label = "TRUE") +
  ggtitle("Aza co-culture, unassigned")
dev.off()

# Labeling clusters for SNU719 and NCC24
new.cluster.ids <- c("SNU719", "1", "2", "NCC24", "4", "5",
                     "6", "7")
names(new.cluster.ids) <- levels(query.dataset)
query.dataset <- RenameIdents(query.dataset, new.cluster.ids)

png("UMAP of labeled Aza co-culture.png")
DimPlot(query.dataset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() +
  ggtitle("Aza co-culture, assigned SNU719 and NCC24 clusters")
dev.off()

# Store Idents (with clusters labeled) as new metadata in coculture.Aza.assigned
coculture.Aza.assigned <- query.dataset
table(Idents(coculture.Aza.assigned))
coculture.Aza.assigned$CellLine <- Idents(coculture.Aza.assigned)

# Store Idents (with clusters labeled) as new metadata in coculture.Aza.assigned
coculture.Aza.assigned <- query.dataset
table(Idents(coculture.Aza.assigned))
coculture.Aza.assigned$CellLine <- Idents(coculture.Aza.assigned)

Idents(coculture.Aza.assigned)<-coculture.Aza.assigned$CellLine
coculture.Aza.assigned$CellLine.Treatment<- paste(Idents(coculture.Aza.assigned), coculture.Aza.assigned$Treatment, sep = "_")
# CellLine.Culture
Idents(coculture.Aza.assigned)<-coculture.Aza.assigned$CellLine
coculture.Aza.assigned$CellLine.Culture <- paste(Idents(coculture.Aza.assigned), coculture.Aza.assigned$Culture, sep = "_")
# CellLine.Treatment.Culture
Idents(coculture.Aza.assigned)<-coculture.Aza.assigned$CellLine.Treatment
coculture.Aza.assigned$CellLine.Treatment.Culture <- paste(Idents(coculture.Aza.assigned), coculture.Aza.assigned$Culture, sep = "_")



save(coculture.Aza.assigned, file="coculture-Aza-assigned.RData")

png("UMAP of assigned Aza-treated SNU719 and NCC24 in co-culture other clusters grey.png")
DimPlot(coculture.Aza.assigned, reduction = "umap", group.by = "CellLine",
        cols = c('SNU719' = 'tomato2',"1" = "gray70", "2" = "gray70", 'NCC24' = 'skyblue2',"4" = "gray70", "5" = "gray70","6" = "gray70", "7" = "gray70")) +
  ggtitle("Aza co-culture, assigned SNU719 and NCC24 clusters")
dev.off()

plot <- DimPlot(coculture.Aza.assigned, reduction = "umap", group.by = "CellLine",
                order = c("SNU719", "1", "2", "NCC24", "4", "5", "6", "7"),
                cols = c('SNU719' = 'tomato2',"1" = "gray70", "2" = "gray70", 'NCC24' = 'gray70',"4" = "gray70", "5" = "gray70","6" = "gray70", "7" = "gray70")) +
  ggtitle("Aza co-culture, assigned SNU719 and NCC24 clusters")

# select.cells <- CellSelector(plot = plot)
# Idents(coculture.Aza.assigned, cells = select.cells) <- "SNU719_filtered"
# coculture.Aza.assigned$CellLine <- Idents(coculture.Aza.assigned)
# save(coculture.Aza.assigned, file="coculture-Aza-assigned-filtered.RData")
