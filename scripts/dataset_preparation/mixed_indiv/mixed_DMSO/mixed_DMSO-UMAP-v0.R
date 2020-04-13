#!/usr/bin/env Rscript

# README ------------------------------------------------------------------

# Script version: 0
# Goal: Plot UMAP of mixed-mixed_DMSO


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

# Load SCtransform-ed and RunPCA-ed object from script mixed_DMSO-v0.R
load(file="/Users/serena/scrnaseq/results/mixed_indiv/mixed_DMSO/v0/SCT-mixed_DMSO-v0.RData")
setwd("/Users/serena/scrnaseq/results/mixed_indiv/mixed_DMSO/v0")

# To run UMAP --------
mixed.DMSO.seuratobject <- RunUMAP(mixed.DMSO.seuratobject, reduction = "pca",dims = 1:9)

mixed.DMSO.seuratobject <- FindNeighbors(mixed.DMSO.seuratobject, dims = 1:15)

mixed.DMSO.seuratobject <- FindClusters(mixed.DMSO.seuratobject,
                                resolution = c(0.4))

## Assign identity of clusters
Idents(object = mixed.DMSO.seuratobject) <- "SCT_snn_res.0.6"

# UMAP Plotting
DimPlot(mixed.DMSO.seuratobject, reduction = "umap", label = "TRUE") + ggtitle("Co-Cultured DMSO-treated Cells") +
  theme(axis.title = element_text(size=12),
        axis.text.x = element_text(size=12,angle = 0, hjust = 0.5),
        axis.text.y = element_text(size=12),
        plot.title=element_text(size=12,face="bold", hjust = 0.5))
beep(20)


# Finding differentially expressed features (cluster biomarkers) ----------------

# Exploring markers
FeaturePlot(mixed.DMSO.seuratobject, features = c("TMSB10", "TP53", "RPS11", "CXCL2", "RPL8", "KRT19", "LYZ", "IGFBP3",
                                                 "CD8A"))
VlnPlot(mixed.DMSO.seuratobject, features = c("IGFBP3"), slot = "counts", assay = "SCT", log = TRUE)

# Finding DE Genes for every cluster compared to all remaining cells. Report only positive ones.
mixed.DMSO.seuratobject.markers <- FindAllMarkers(mixed.DMSO.seuratobject, min.pct = 0.25, logfc.threshold = 0.25)
mixed.DMSO.seuratobject.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# Combine markers with gene descriptions
# Report top10 DE genes of each cluster compared to all remainging cells.
mixed.DMSO.seuratobject.markers_ann <- mixed.DMSO.seuratobject.markers  %>%
  rownames_to_column(var="genes") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name")) %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_logFC)
write.csv(mixed.DMSO.seuratobject.markers_ann,"top10allmarkersann.csv")

save(mixed.DMSO.seuratobject, file="/Users/serena/scrnaseq/results/mixed_indiv/mixed_DMSO/v0/UMAP-Dims9-Res06-mixed-DMSO.Rdata")
save(mixed.DMSO.seuratobject.markers, file="/Users/serena/scrnaseq/results/mixed_indiv/mixed_DMSO/v0/DEAllclusters-Dims9-Res06-mixed-DMSO.Rdata")
save(mixed.DMSO.seuratobject.markers_ann, file="/Users/serena/scrnaseq/results/mixed_indiv/mixed_DMSO/v0/ANNDEAllclusters-Dims9-Res06-mixed-DMSO.Rdata")
beep(2)

# Plotting volcano plots for each individual cluster----------------------------------------
#Load annotations
annotations <- read.csv("/users/serena/scrnaseq/data/annotation.csv")

### Find all markers of Cluster 1
cluster1.markers <- FindMarkers(mixed.DMSO.seuratobject, ident.1 = 1, min.pct = 0.25)


# Combine markers with gene descriptions
cluster1.markers_ann <- cluster1.markers %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
write.csv(cluster1.markers_ann,"cluster1markersann.csv")

## Sort by ordered padj
cluster1.markers_ann_ordered <- cluster1.markers_ann[order(cluster1.markers_ann$p_val_adj), ]

## Create a column to indicate which genes to label
cluster1.markers_ann_ordered$genelabels <- ""
cluster1.markers_ann_ordered$genelabels[1:10] <- T

# Create column of -log10(p_val_adj)
`-log10(p_val_adj)` <- -log10(cluster1.markers_ann_ordered$p_val_adj)
cbind(cluster1.markers_ann_ordered,`-log10(p_val_adj)`)

# Set threshold. Inspect .csv file first.
threshold <- cluster1.markers_ann_ordered$`p_val_adj` < 10e-25

## Add logical vector as a column (threshold) to the data frame
cluster1.markers_ann_ordered$threshold <- threshold

# Plot labeled Volcano Plot
ggplot(cluster1.markers_ann_ordered, aes(x=`avg_logFC`/log(2), y=(-log10(`p_val_adj`)), color=ifelse(cluster1.markers_ann_ordered$threshold == T, "grey","red"))) +
  geom_point() +
  geom_text_repel(aes(x = avg_logFC/log(2), y = -log10(p_val_adj), label = ifelse(cluster1.markers_ann_ordered$genelabels == T, cluster1.markers_ann_ordered$gene,""))) +
  ggtitle("Volcano Plot SNU719-DMSO-april2017 Cluster 1") +
  xlab("log2 FC") +
  ylab("-log10 adjusted p-value") +
  # scale_y_continuous(limits = c(0,50)) +
  scale_x_continuous(limits=c(-2,2)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

# # find all markers distinguishing cluster 4 from clusters 0 and 2
# cluster4.markers <- FindMarkers(mixed.DMSO.seuratobject, ident.1 = 4, ident.2 = c(0, 3), min.pct = 0.25)

# # find markers for every cluster compared to all remaining cells, report only the positive ones
# SNU719.DMSO.april2017.markers <- FindAllMarkers(mixed.DMSO.seuratobject, min.pct = 0.25, logfc.threshold = 0.25)
# SNU719.DMSO.april2017.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


# Generate heatmap
# Plotting the top 20 markers (or all markers if less than 20) for each cluster.
# First: find markers for every cluster compared to all remaining cells, report only the positive ones
# mixed.DMSO.seuratobject.markers <- FindAllMarkers(mixed.DMSO.seuratobject, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# #change n to chaneg number of markers reported.
# mixed.DMSO.seuratobject.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
# top10 <- mixed.DMSO.seuratobject.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# png("HeatmapDEgenesofeachCluster.png")
# DoHeatmap(mixed.DMSO.seuratobject, features = top10$gene) + NoLegend()
# dev.off()