#!/usr/bin/env Rscript

# README ------------------------------------------------------------------
# Assign some cells in DMSO-treated co-culture as SNU719 or NCC24,
# using the monoculture DMSO-treated SNU719 and monoculture DMSO-treated NCC24 cells as a reference.

# Set-up packages, working directory, datasets -----------------------------------------------------------
## Load Packages
library(Seurat)
library(sctransform)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
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
setwd("/handover/results/cell line assignment/DMSO")

## Load datasets
# Load DMSO-treated SNU719 in monoculture and DMSO-treated NCC24 in monoculture for building reference dataset
load("/handover/results/dataset preparation/mono/april2017/SNU719-DMSO-april2017/SCT-SNU719-DMSO-april2017.RData")
load("/handover/results/dataset preparation/mono/july2017/SNU719-DMSO-july2017/SCT-SNU719-DMSO-july2017.RData")
load("/handover/results/dataset preparation/mono/sept2017/NCC24-DMSO-sept2017/SCT-NCC24-DMSO-sept2017.RData")
# Load DMSO-treated co-culture
load("/handover/results/dataset preparation/co/co-DMSO/SCT-co_DMSO.RData")

# Forming reference dataset -----------------------------------------------------------
# Merge SNU719 and NCC24  to form a reference (object: monoculture.DMSO.merged)
monoculture.DMSO.merged <- merge(x = NCC24.DMSO.sept2017.seuratobject,
                          y = c(SNU719.DMSO.april2017.seuratobject,SNU719.DMSO.july2017.seuratobject),
                          add.cell.id = c("NCC24.DMSO.sept2017", "SNU719.DMSO.april2017", "SNU719.DMSO.july2017"),
                          project = "merge_ref_DMSO",
                          merge.data = TRUE)
# Run SCTransform
monoculture.DMSO.merged <-SCTransform(monoculture.DMSO.merged,
                               return.only.var.genes = F,
                               vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
# Run PCA and plot Elbow Plot
monoculture.DMSO.merged <- RunPCA(monoculture.DMSO.merged, verbose = FALSE)
png("Elbow Plot of merged DMSO-treated SNU719 and NCC24 in monoculture.png")
ElbowPlot(monoculture.DMSO.merged) + ggtitle("Merged DMSO-treated SNU719 and NCC24 in monoculture")
dev.off()

png("PCA Plot of merged DMSO-treated SNU719 and NCC24 in monoculture.png")
PCAPlot(monoculture.DMSO.merged, group.by = "sample") +
  ggtitle("Merged DMSO-treated SNU719 and NCC24 in monoculture")
dev.off()

monoculture.DMSO.merged <- RunUMAP(monoculture.DMSO.merged, reduction = "pca",dims = 1:15)
monoculture.DMSO.merged <- FindNeighbors(monoculture.DMSO.merged, dims = 1:15)
monoculture.DMSO.merged <- FindClusters(monoculture.DMSO.merged, resolution = c(0.05))

# Observe clustering with UMAP at resolution of 0.05
Idents(object = monoculture.DMSO.merged) <- "SCT_snn_res.0.05"
png("UMAP of merged DMSO-treated SNU719 and NCC24 in monoculture.png")
DimPlot(monoculture.DMSO.merged, reduction = "umap", group.by = "CellLine",
        cols = c('SNU719' = 'tomato2','NCC24' = 'skyblue2')) +
  ggtitle("Merged DMSO-treated SNU719 and NCC24 in monoculture")
dev.off()

# Save reference dataset
save(monoculture.DMSO.merged, file="monoculture_DMSO_merged.RData")



# Cell line assignment in co-culture -----------------------------------------------------------
reference.dataset <- monoculture.DMSO.merged
query.dataset <- coculture.DMSO

DefaultAssay(query.dataset) <- "SCT"

# Set reduction to CCA for FindTransferAnchors
query.anchors <- FindTransferAnchors(reference = reference.dataset, query = query.dataset, dims = 1:10,
                                     reduction = 'cca')
predictions <- TransferData(anchorset = query.anchors, refdata = reference.dataset$CellLine, dims = 1:10,
                            weight.reduction = 'cca')
query.dataset <- AddMetaData(query.dataset, metadata = predictions)

# Save predicted.id results to a csv file
write.csv(query.dataset@meta.data[["predicted.id"]], file = "Predicted Cell Line Assignment DMSO.csv")
save(predictions, file = "predictions.Rdata")
# Plotting prediction scores of SNU719 against NCC24 to observe the rs------------------------------------------------------------
# Remarks: They are inversely linearly related.
plot(predictions$prediction.score.NCC24,predictions$prediction.score.SNU719)
ggplot(predictions, aes(x=prediction.score.NCC24, y=prediction.score.SNU719)) + geom_point(size = 0.1, shape = 23)
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
query.dataset.sc1 <- FindClusters(query.dataset.sc1, resolution = c(0.4), dims = 1:15)
Idents(object = query.dataset.sc1) <- "SCT_snn_res.0.4"

png("UMAP of top 1500 cells likely to be SNU719 and NCC24 for DMSO co-culture.png")
DimPlot(query.dataset.sc1, reduction = "umap", group.by = "predicted.id.sc1", label = "FALSE",
        cols = c('SNU719' = 'tomato2','NCC24' = 'skyblue2', 'Unassigned' = 'gray70')) +
  ggtitle("Top 1500 cells likely to be SNU719 and NCC24 for DMSO co-culture")
dev.off()

# png("UMAP of top 1500 cells likely to be SNU719 and NCC24 for DMSO co-culture on same axis as monoculture.png")
# DimPlot(query.dataset.sc1, reduction = "umap", group.by = "predicted.id.sc1", label = "FALSE",
#         cols = c('SNU719' = 'tomato2','NCC24' = 'skyblue2', 'Unassigned' = 'gray70')) +
#   xlim(-8,4) +
#   ylim(-2.5,5) +
#   ggtitle("Top 1500 cells likely to be SNU719 and NCC24 for DMSO co-culture")
# dev.off()

# Labeling cells in co-culture as SNU719 and NCC24------------------------------------------------------------------
query.dataset <- RunUMAP(query.dataset, reduction = "pca", dims = 1:15)
query.dataset <- FindNeighbors(query.dataset, dims = 1:15)
query.dataset <- FindClusters(query.dataset, resolution = c(0.4), dims = 1:15)
# Observe clustering at resolution of 0.4
Idents(object = query.dataset) <- "SCT_snn_res.0.4"

png("UMAP of cells likely to be SNU719 and NCC24 for DMSO co-culture.png")
DimPlot(query.dataset, reduction = "umap", group.by = "predicted.id", label = "TRUE",
        cols = c('SNU719' = 'tomato2','NCC24' = 'skyblue2')) +
  ggtitle("If all cells are SNU719 or NCC24 for DMSO co-culture")
dev.off()

png("UMAP of unlabeled DMSO co-culture.png")
DimPlot(query.dataset, reduction = "umap", label = "TRUE") +
  ggtitle("DMSO co-culture, unassigned")
dev.off()

# Labeling clusters for SNU719 and NCC24
new.cluster.ids <- c("SNU719", "1", "2", "NCC24", "4", "5",
                     "6", "7")
names(new.cluster.ids) <- levels(query.dataset)
query.dataset <- RenameIdents(query.dataset, new.cluster.ids)

png("UMAP of labeled DMSO co-culture.png")
DimPlot(query.dataset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() +
  ggtitle("DMSO co-culture, assigned SNU719 and NCC24 clusters")
dev.off()

# Store Idents (with clusters labeled) as new metadata in coculture.DMSO.assigned
coculture.DMSO.assigned <- query.dataset
table(Idents(coculture.DMSO.assigned))
coculture.DMSO.assigned$CellLine <- Idents(coculture.DMSO.assigned)

save(coculture.DMSO.assigned, file="coculture-DMSO-assigned.RData")


png("UMAP of assigned DMSO-treated SNU719 and NCC24 in co-culture other clusters grey.png")
DimPlot(coculture.DMSO.assigned, reduction = "umap", group.by = "CellLine",
        cols = c('SNU719' = 'tomato2',"1" = "gray70", "2" = "gray70", 'NCC24' = 'skyblue2',"4" = "gray70", "5" = "gray70","6" = "gray70", "7" = "gray70")) +
  ggtitle("DMSO co-culture, assigned SNU719 and NCC24 clusters")
dev.off()


# plot <- DimPlot(coculture.DMSO.assigned, reduction = "umap", group.by = "CellLine",
#                 order = c("SNU719", "1", "2", "NCC24", "4", "5", "6", "7"),
#                 cols = c('SNU719' = 'tomato2',"1" = "gray70", "2" = "gray70", 'NCC24' = 'gray70',"4" = "gray70", "5" = "gray70","6" = "gray70", "7" = "gray70")) +
#   ggtitle("DMSO co-culture, assigned SNU719 and NCC24 clusters")

# select.cells <- CellSelector(plot = plot)
# Idents(coculture.DMSO.assigned, cells = select.cells) <- "SNU719_filtered"
# coculture.DMSO.assigned$CellLine <- Idents(coculture.DMSO.assigned)
# save(coculture.DMSO.assigned, file="/scrnaseq/data/coculture-DMSO-assigned-filtered.RData")
