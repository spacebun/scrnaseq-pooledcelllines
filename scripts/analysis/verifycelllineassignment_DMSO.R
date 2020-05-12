#!/usr/bin/env Rscript

# README ------------------------------------------------------------------
#  Verifying accuracy of cell line assignment of SNU719 and NCC24 in DMSO-treated co-culture .

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
setwd("/scrnaseq/results/analysis/verify cell line assignment/DMSO")

## Load datasets
# Load SNU719 and NCC24 DMSO-treated cells in monoculture
load("/scrnaseq/results/cell line assignment/DMSO/monoculture_DMSO_merged.RData")
# Load co-culture with assigned SNU719 and DMSO cells
load("/scrnaseq/results/cell line assignment/DMSO/coculture-DMSO-assigned.RData")

# Transfer object to new object for analysis
snu.ncc <- monoculture.DMSO.merged

# Subset SNU719 and NCC24 labeled cells in coculture
Idents(coculture.DMSO.assigned) <- coculture.DMSO.assigned$CellLine
coculture.DMSO.cells <- subset(coculture.DMSO.assigned, idents=c("NCC24","SNU719"))
DefaultAssay(coculture.DMSO.cells) <- "RNA"
coculture.DMSO.cells <- NormalizeData(coculture.DMSO.cells)
coculture.DMSO.cells <- ScaleData(coculture.DMSO.cells)
DefaultAssay(coculture.DMSO.cells) <- "SCT"

# DE between both cell lines-----------------------------------------------------------
# Switch current ident to column holding CellLine metadata
Idents(snu.ncc)<-snu.ncc$CellLine
annotations <- read.csv("/scrnaseq/data/annotation.csv")

# FindAllMarkers
CellLine.markers.FindAll <- FindAllMarkers(snu.ncc, idents= 'CellLine', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CellLine.markers.FindAll %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# Combine markers with gene descriptions
CellLine.markers.FindAll_ann <- CellLine.markers.FindAll%>%
  rownames_to_column(var="genes") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

write.csv(CellLine.markers.FindAll_ann,"Markers for SNU719 vs NCC24 in monoculture.csv")

# FindAllMarkers: Subset top 10 most DE markers for each Cell Line
CellLine.markers.FindAll.heatmap <- CellLine.markers.FindAll_ann %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_logFC) %>%
  # rownames_to_column(var="genes") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# Plot heatmap of top 10 DE genes (from FindAllMarkers) against Cell Lines
Idents(snu.ncc)<-snu.ncc$CellLine
levels(snu.ncc)
levels(snu.ncc) <- c("SNU719","NCC24")
levels(snu.ncc)

DoHeatmap(snu.ncc, features = CellLine.markers.FindAll.heatmap$gene, group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
  scale_fill_distiller(palette = "RdBu") + ggtitle("DMSO-treated SNU719 and DMSO-treated NCC24 cells in monoculture")


DoHeatmap(snu.ncc, features = CellLine.markers.FindAll.heatmap$gene, group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
  scale_fill_distiller(palette = "RdBu") + ggtitle("DMSO-treated SNU719 and DMSO-treated NCC24 cells in monoculture")


# Plot the same DEG on the mixed.DMSO.cells
DoHeatmap(coculture.DMSO.cells, features = CellLine.markers.FindAll.heatmap$gene, group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
  scale_fill_distiller(palette = "RdBu") + ggtitle("DMSO-treated SNU719 and DMSO-treated NCC24 cells in co-culture")
