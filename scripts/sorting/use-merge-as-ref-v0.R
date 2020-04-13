#!/usr/bin/env Rscript

# README ------------------------------------------------------------------
# Script version: 0
# Goal: To use merged DMSO lines as reference dataset and query dataset

# Outputs from this script:

# Integrated SNU719-DMSO objects (SNU.DMSO.Integrated)
# 

# Merged SNU and NCC object (pure.DMSO.merged/reference.dataset)
# 

# Labeled mixed_DMSO query dataset (mixed_DMSO_labeled)
# Dims 15 Res 04
#  NCC24 SNU719 
#  1811    3777
#  C3       C0
# 

# Steps Done:
# Formed reference dataset 
# Labeled all cells in query dataset as either SNU719 or NCC24 (query.dataset)
# Labeled top 1500 cells with highest prediction scores for each cell line in the query dataset (query.dataset.sc1)
# Forced Cluster Query dataset with 0.4 (7 clusters) or 0.8 (6 clusters)
# Performed FindAllMarkers

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
setwd("/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0")

## Load the SCT Datasets
load("/Users/serena/scrnaseq/results/pure_indiv/april2017/SNU719-DMSO-april2017/v1/SCT-SNU719-DMSO-april2017-v1.RData")
load("/Users/serena/scrnaseq/results/pure_indiv/july2017/SNU719-DMSO-july2017/v1/SCT-SNU719-DMSO-july2017-v1.RData")
load("/Users/serena/scrnaseq/results/pure_indiv/sept2017/NCC24-DMSO-sept2017/v1/SCT-NCC24-DMSO-sept2017-v1.RData")
load("/Users/serena/scrnaseq/results/mixed_indiv/mixed_DMSO/v0/SCT-mixed_DMSO-v0.RData")

# Forming reference dataset -----------------------------------------------------------
# Part 1 of 2: Integrating SNU719-DMSO Datasets
# This may not be necessary as merging in Part 2 does not merge with the integrated counts. Merge with raw counts?
SNU.DMSO.list <- list( SNU719.DMSO.april2017.seuratobject,SNU719.DMSO.july2017.seuratobject)
SNU.DMSO.features <- SelectIntegrationFeatures(object.list = SNU.DMSO.list, nfeatures = 3000)
SNU.DMSO.list <- PrepSCTIntegration(object.list = SNU.DMSO.list, anchor.features = SNU.DMSO.features)
SNU.DMSO.anchors <- FindIntegrationAnchors(object.list = SNU.DMSO.list, normalization.method = "SCT", anchor.features = SNU.DMSO.features)
SNU.DMSO.integrated <- IntegrateData(anchorset = SNU.DMSO.anchors, normalization.method = "SCT")
DefaultAssay(SNU.DMSO.integrated) <- "integrated"

save(SNU.DMSO.integrated, file="/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/DMSO Reference/integrate_SNU_DMSO.RData")
# load("/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/integrate_SNU_DMSO.RData")

# Part 2 of 2: Merging SNU719 and NCC24 datasets to form a reference
pure.DMSO.merged <- merge(x = NCC24.DMSO.sept2017.seuratobject,
                y = c(SNU.DMSO.integrated),
                       add.cell.id = c("NCC24.DMSO.sept2017","SNU719.DMSO.integrated"),
                     project = "merge_ref_DMSO",
                merge.data = TRUE)
pure.DMSO.merged <-SCTransform(pure.DMSO.merged,
                          return.only.var.genes = F,
                          vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
pure.DMSO.merged <- RunPCA(pure.DMSO.merged, verbose = FALSE)
# Use Elbow Plot to assess dimensionality of dataset for RunUMAP. 
# Remarks: Elbow around 15
ElbowPlot(pure.DMSO.merged)
# Observe PCAPlot of object by sample behaviour
PCAPlot(pure.DMSO.merged,
        group.by = "sample")
pure.DMSO.merged <- RunUMAP(pure.DMSO.merged, reduction = "pca",dims = 1:15)
pure.DMSO.merged <- FindNeighbors(pure.DMSO.merged, dims = 1:15)
pure.DMSO.merged <- FindClusters(pure.DMSO.merged, resolution = c(0.05,0.07,0.1,0.6,1.0,1.6,2.0))
# Observe clustering at resolution of 0.05
Idents(object = pure.DMSO.merged) <- "SCT_snn_res.0.05"
# Plot UMAP by clusters
DimPlot(pure.DMSO.merged, reduction = "umap",
        cols = c('0' = 'tomato2', '1' = 'skyblue2'))
# Plot UMAP by defined Cell Line
DimPlot(pure.DMSO.merged, reduction = "umap", group.by = "CellLine",
        cols = c('SNU719' = 'tomato2','NCC24' = 'skyblue2'))
DimPlot(pure.DMSO.merged, reduction = "umap", group.by = "orig.ident")

# Obtain number of cells in each cluster
table(pure.DMSO.merged@meta.data$SCT_snn_res.0.05)
# Remarks: Cluster 0 (red) corresponds to SNU719
# Remarks: Cluster 1 (blue) corresponds to NCC24

# Save reference dataset
save(pure.DMSO.merged, file="/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/merged_SNUintegrated_NCC_DMSO.RData")
# load("/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/merged_SNUintegrated_NCC_DMSO.RData")

# Using the merged datasets as a reference to query cells in mixed dataset----------------------
reference.dataset <- pure.DMSO.merged
query.dataset <- mixed.DMSO.seuratobject

# Investigate Elbow Plot of query dataset. 
# Remarks: Elbow around 10
ElbowPlot(query.dataset)

# Perform querying
# Remarks: Chosen reduction = 'cca' due to small number of anchors found when using PCA
DefaultAssay(query.dataset) <- "SCT"
query.anchors <- FindTransferAnchors(reference = reference.dataset, query = query.dataset, dims = 1:10,
    reduction = 'cca')
predictions <- TransferData(anchorset = query.anchors, refdata = reference.dataset$CellLine, dims = 1:10,
    weight.reduction = 'cca')
query.dataset <- AddMetaData(query.dataset, metadata = predictions)

# Save predicted.id results to a csv file
write.csv(query.dataset@meta.data[["predicted.id"]], file = "Predicted Cell Lines Dims 10 with CCA.csv")
save(query.dataset, file="/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/query_dataset_ref_DMSO_usingCCA_v0.RData")

#### Visualize query dataset on a UMAP, separated by Predicted Cell Line---------------------------------------------------
# Note Chosen dims is 30 -> why? should be 15
query.dataset <- RunUMAP(query.dataset, reduction = "pca", dims = 1:30)
query.dataset <- FindNeighbors(query.dataset, dims = 1:15)
query.dataset <- FindClusters(query.dataset, resolution = c(0.4,0.6,1.0), dims = 1:15)
# Observe clustering at resolution of 0.4
Idents(object = query.dataset) <- "SCT_snn_res.0.4"

DimPlot(query.dataset, reduction = "umap", label = "TRUE")

DimPlot(query.dataset, reduction = "umap", group.by = "predicted.id", label = "TRUE",
        cols = c('SNU719' = 'tomato2','NCC24' = 'skyblue2'))

# Number of cells in each cluster
table(query.dataset@meta.data$SCT_snn_res.0.4)

# Remarks: Cluster 0 == SNU719 and Cluster 3 == NCC24 in query dataset
# Parameters: Dims at 30, Resolution at 0.4

#### Reorganised predicted.id by prediction score, then plot UMAP again.---------------------------------------------------

### Label Top 1500 cells with highest prediction score for each cell line
# Transfer predictions and query dataset to new variables
predictions.sc1 <- predictions
query.dataset.sc1 <- query.dataset
## Create a column to indicate which cells to label
predictions.sc1$predicted.id.sc1 <- "Unassigned"
## Sort by ordered SNU719 prediction scores
predictions.sc1 <-predictions.sc1[order(-predictions.sc1$prediction.score.SNU719), ]
## Label top 1500 cells as SNU719
predictions.sc1$predicted.id.sc1[1:1500] <- "SNU719"
## Sort by ordered NCC24 prediction scores
predictions.sc1 <-predictions.sc1[order(-predictions.sc1$prediction.score.NCC24), ]
## Label top 1500 cells as NCC24
predictions.sc1$predicted.id.sc1[1:1500] <- "NCC24"
## Add metadata back to query.dataset.sc1
query.dataset.sc1 <- AddMetaData(query.dataset.sc1, metadata = predictions.sc1)
## Visualize with UMAP
query.dataset.sc1 <- RunUMAP(query.dataset.sc1, reduction = "pca", dims = 1:30)
query.dataset.sc1 <- FindNeighbors(query.dataset.sc1,
                                   dims = 1:15)
query.dataset.sc1 <- FindClusters(query.dataset.sc1,
                                  resolution = c(0.4,0.6,1.0),
                                  dims = 1:15)

Idents(object = query.dataset.sc1) <- "SCT_snn_res.0.4"

DimPlot(query.dataset.sc1, reduction = "umap", label = "TRUE")

DimPlot(query.dataset.sc1, reduction = "umap", group.by = "predicted.id.sc1", label = "TRUE",
        cols = c('SNU719' = 'tomato2','NCC24' = 'skyblue2', 'Unassigned' = 'gray70'))

save(query.dataset, file="/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/query_dataset_ref_DMSO_usingCCA_v0.RData")


# Find DE markers that differentiate clusters NCC24 and SNU719 in query dataset --------------
# These markers should correspond to the DE markers found between mono datasets SNU719 and NCC24.
# Cross check these markers with DE-mono-v0.R DE genes of mono SNU v NCC.
#Load annotations
annotations <- read.csv("/users/serena/scrnaseq/data/annotation.csv")

# All Cluster 0 (SNU719) Markers
cluster0.markers <- FindMarkers(query.dataset.sc1, ident.1 = 0, min.pct = 0.25)
# Combine markers with gene descriptions
cluster0.markers_ann <- cluster0.markers %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
  # group_by(cluster)
  # top_n(n = 10, wt = avg_logFC)
write.csv(cluster0.markers_ann,"query_cluster0_SNU719_markers_ann.csv")

# Markers distinguishing Cluster 0 (SNU719) from Cluster 3 (NCC24)
cluster0.cluster3.markers <- FindMarkers(query.dataset.sc1, ident.1 = 0, ident.2 = 3, min.pct = 0.25)
# Combine markers with gene descriptions
cluster0.cluster3.markers_ann <- cluster0.cluster3.markers %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
  # group_by(cluster) %>%
  # top_n(n = 10, wt = avg_logFC)
write.csv(cluster0.cluster3.markers_ann,"query_cluster0_v_cluster3_markers_ann.csv")

# Labeling cells in mixed dataset as SNU719 ------------------------------------------------------------------
Idents(object = query.dataset) <- "SCT_snn_res.0.4"

new.cluster.ids <- c("SNU719", "1", "2", "NCC24", "4", "5", 
                     "6", "7")
names(new.cluster.ids) <- levels(query.dataset)
query.dataset <- RenameIdents(query.dataset, new.cluster.ids)
DimPlot(query.dataset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

mixed_DMSO_labeled <- query.dataset
# Store Idents (with clusters labeled) as new metadata
table(Idents(mixed_DMSO_labeled))
mixed_DMSO_labeled$CellLine <- Idents(mixed_DMSO_labeled)

save(mixed_DMSO_labeled, file="/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/mixed_DMSO_labeled_v0.RData")


