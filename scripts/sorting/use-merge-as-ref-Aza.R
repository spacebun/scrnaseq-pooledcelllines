Aza#!/usr/bin/env Rscript

# README ------------------------------------------------------------------

# Script version: 0
# Goal: To identify

# Output from this script:
# Integrated SNU719-Aza objects
# /Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/integrate_SNU_Aza.RData

# Merged SNU and NCC object (reference dataset)
# /Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/merged_SNUintegrated_NCC_Aza.RData

# Labeled mixed_Aza query dataset (all labeled either SNU719 or NCC24)
#  NCC24 SNU719
# 2442   1059


# Labeled mixed_DMSO query dataset 
# /Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/mixed_Aza_labeled_v0.RData

# UMAP

# Steps Done:
# Formed reference dataset (reference.dataset or pure.DMSO.merged)
# Labeled all cells in query dataset as either SNU719 or NCC24 (query.dataset)
# Labeled top 1500 cells with highest prediction scores for each cell line in the query dataset (query.dataset.sc1)
# Forced Cluster Query dataset with 0.4 (6 clusters) 
# 

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

## Load the Datasets
load("/Users/serena/scrnaseq/results/mixed_indiv/mixed_Aza/v0/SCT-mixed_Aza-v0.RData")
load("/Users/serena/scrnaseq/results/mixed_indiv/mixed_DMSO/v0/SCT-mixed_DMSO-v0.RData")

load("/Users/serena/scrnaseq/results/pure_indiv/april2017/SNU719-Aza-april2017/v1/SCT-SNU719-Aza-april2017-v1.RData")
load("/Users/serena/scrnaseq/results/pure_indiv/april2017/SNU719-DMSO-april2017/v1/SCT-SNU719-DMSO-april2017-v1.RData")

load("/Users/serena/scrnaseq/results/pure_indiv/july2017/SNU719-DMSO-july2017/v1/SCT-SNU719-DMSO-july2017-v1.RData")
# load("/Users/serena/scrnaseq/results/pure_indiv/july2017/HFE145-july2017/v1/SCT-HFE145-july2017-v1.RData")

load("/Users/serena/scrnaseq/results/pure_indiv/sept2017/SNU719-Aza-sept2017/v1/SCT-SNU719-Aza-sept2017-v1.RData")
load("/Users/serena/scrnaseq/results/pure_indiv/sept2017/NCC24-DMSO-sept2017/v1/SCT-NCC24-DMSO-sept2017-v1.RData")
# load("/Users/serena/scrnaseq/results/pure_indiv/sept2017/YCC10-Aza-sept2017/v1/SCT-YCC10-Aza-sept2017-v1.RData")
# load("/Users/serena/scrnaseq/results/pure_indiv/sept2017/YCC10-DMSO-sept2017/v1/SCT-YCC10-DMSO-sept2017-v1.RData")

# Integrate SNU719 Datasets to create one Seurat Object -----------------------------------------------------------
SNU.Aza.list <- list( SNU719.Aza.april2017.seuratobject,SNU719.Aza.sept2017.seuratobject)
SNU.Aza.features <- SelectIntegrationFeatures(object.list = SNU.Aza.list, nfeatures = 3000)
SNU.Aza.list <- PrepSCTIntegration(object.list = SNU.Aza.list, anchor.features = SNU.Aza.features)
SNU.Aza.anchors <- FindIntegrationAnchors(object.list = SNU.Aza.list, normalization.method = "SCT", anchor.features = SNU.Aza.features)
beep(1)
SNU.Aza.integrated <- IntegrateData(anchorset = SNU.Aza.anchors, normalization.method = "SCT")
beep(8)
DefaultAssay(SNU.Aza.integrated) <- "integrated"
# Save integrated seurat object
save(SNU.Aza.integrated, file="/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/integrate_SNU_Aza.RData")
#load("/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/integrate_SNU_Aza.RData")
beep(2)

# Merge SNU & NCC,plot UMAP to obtain 2 clusters that correspond to each cell line-----------------------------------------------------------
pure.Aza.merged <- merge(x = NCC24.DMSO.sept2017.seuratobject,
                y = c(SNU.Aza.integrated),
                       add.cell.id = c("NCC24.DMSO.sept2017","SNU719.Aza.integrated"),
                     project = "merge_ref_Aza")

pure.Aza.merged <-SCTransform(pure.Aza.merged,
                          return.only.var.genes = F,
                          vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
beep(2)
pure.Aza.merged <- RunPCA(pure.Aza.merged, verbose = FALSE)

#Assess dimensionality of dataset for RunUMAP
# Remarks: Elbow around 15
ElbowPlot(pure.Aza.merged)

# Observe PCAPlot of object by sample behaviour
PCAPlot(pure.Aza.merged,
        group.by = "sample")

# Adjust dims and resolution according to Elbow Plot and number of cells
# Chosen dims: 15
pure.Aza.merged <- RunUMAP(pure.Aza.merged, reduction = "pca",dims = 1:15)
pure.Aza.merged <- FindNeighbors(pure.Aza.merged, dims = 1:15)
pure.Aza.merged <- FindClusters(pure.Aza.merged,
                                resolution = c(0.05,0.07,0.1,0.2,0.5))
# After adjustment with FindClusters and DimPlot, chosen resn = 0.1
Idents(object = pure.Aza.merged) <- "SCT_snn_res.0.1"

DimPlot(pure.Aza.merged, reduction = "umap",
        cols = c('0' = 'tomato2', '1' = 'skyblue2'))

DimPlot(pure.Aza.merged, reduction = "umap", group.by = "CellLine",
        cols = c('SNU719' = 'tomato2','NCC24' = 'skyblue2'))

DimPlot(pure.Aza.merged, reduction = "umap", group.by = "sample")
# Number of cells in each cluster
table(pure.Aza.merged@meta.data$SCT_snn_res.0.1)

# Remarks: Cluster 0 (red) corresponds to SNU719
# Remarks: Cluster 1 (blue) corresponds to NCC24

# pure.DMSO.merged <- RunTSNE(pure.DMSO.merged)
# DimPlot(pure.DMSO.merged,reduction = "tsne")

# Save merged object of SNU719 & NCC24
save(pure.Aza.merged, file="/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/merged_SNUintegrated_NCC_Aza.RData")
# load("/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/merged_SNUintegrated_NCC_Aza.RData")

# Using the merged datasets as a reference to query mixed cells----------------------
reference.dataset <- pure.Aza.merged
query.dataset <- mixed.Aza.seuratobject

# Investigate Elbow Plot
# Remarks: Elbow around 10
ElbowPlot(reference.dataset)
ElbowPlot(query.dataset)

# Adjust dims to 10 based on ElbowPlots
# Chosen reduction = 'cca' due to small number of anchors found when using PCA
DefaultAssay(query.dataset) <- "SCT"
query.anchors <- FindTransferAnchors(reference = reference.dataset, query = query.dataset,
    dims = 1:10,
    reduction = 'cca')
beep(1)
predictions <- TransferData(anchorset = query.anchors, refdata = reference.dataset$CellLine,
    dims = 1:10,
    weight.reduction = 'cca')
query.dataset <- AddMetaData(query.dataset, metadata = predictions)

# Save predicted.id results to a csv file
write.csv(query.dataset@meta.data[["predicted.id"]], file = "Predicted Cell Lines Dims 10 with CCA for Aza.csv")
save(query.dataset, file="/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/query_dataset_ref_Aza_usingCCA_v0.RData")
beep(1)

#### Visualize query dataset on a UMAP, separated by Predicted Cell Line---------------------------------------------------
# Note Chosen dims is 30 -> why? should be 15
query.dataset <- RunUMAP(query.dataset, reduction = "pca", dims = 1:30)
query.dataset <- FindNeighbors(query.dataset,
                                     dims = 1:15)
query.dataset <- FindClusters(query.dataset,
                                    resolution = c(0.4,0.5,0.6,1.0),
                                    dims = 1:15)
# After adjustment with FindClusters and DimPlot, chosen resn = 0.4
# This gives approx 7 clusters
Idents(object = query.dataset) <- "SCT_snn_res.0.5"

DimPlot(query.dataset, reduction = "umap", label = "TRUE")

DimPlot(query.dataset, reduction = "umap", group.by = "predicted.id", label = "TRUE",
        cols = c('SNU719' = 'tomato2','NCC24' = 'skyblue2'))

# Number of cells in each cluster
table(query.dataset@meta.data$SCT_snn_res.0.6)

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
predictions.sc1$predicted.id.sc1[1:700] <- "NCC24"
## Add metadata back to query.dataset.sc1
query.dataset.sc1 <- AddMetaData(query.dataset.sc1, metadata = predictions.sc1)
## Visualize with UMAP
query.dataset.sc1 <- RunUMAP(query.dataset.sc1, reduction = "pca", dims = 1:30)
query.dataset.sc1 <- FindNeighbors(query.dataset.sc1,
                                   dims = 1:15)
query.dataset.sc1 <- FindClusters(query.dataset.sc1,
                                  resolution = c(0.4,0.5,0.6),
                                  dims = 1:15)

Idents(object = query.dataset.sc1) <- "SCT_snn_res.0.5"

DimPlot(query.dataset.sc1, reduction = "umap", label = "TRUE")

DimPlot(query.dataset.sc1, reduction = "umap", group.by = "predicted.id.sc1", label = "TRUE",
        cols = c('SNU719' = 'tomato2','NCC24' = 'skyblue2', 'Unassigned' = 'gray70'))

save(query.dataset, file="/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/query_dataset_ref_Aza_usingCCA_v0.RData")


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
Idents(object = query.dataset) <- "SCT_snn_res.0.5"

new.cluster.ids <- c("SNU719", "1", "2", "NCC24", "4", "5",
                     "6", "7")
names(new.cluster.ids) <- levels(query.dataset)
query.dataset <- RenameIdents(query.dataset, new.cluster.ids)
DimPlot(query.dataset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

mixed_Aza_labeled <- query.dataset

# Store Idents (with clusters labeled) as new metadata
table(Idents(mixed_Aza_labeled))
mixed_Aza_labeled$CellLine <- Idents(mixed_Aza_labeled)

save(mixed_Aza_labeled, file="/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/mixed_Aza_labeled_v0.RData")
