#!/usr/bin/env Rscript

# README ------------------------------------------------------------------
#  Verifying accuracy of cell line assignment of SNU719 and NCC24 in Aza-treated co-culture .

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
setwd("/scrnaseq/results/analysis/verify cell line assignment/Aza")

## Load datasets
# Load SNU719 and NCC24 Aza-treated cells in monoculture
load("/scrnaseq/results/cell line assignment/Aza/monoculture_Aza_merged.RData")
# Load co-culture with assigned SNU719 and Aza cells
load("/scrnaseq/results/cell line assignment/Aza/coculture-Aza-assigned.RData")

# Transfer object to new object for analysis
snu.ncc <- monoculture.Aza.merged

# Subset SNU719 and NCC24 labeled cells in coculture
Idents(coculture.Aza.assigned) <- coculture.Aza.assigned$CellLine
coculture.Aza.cells <- subset(coculture.Aza.assigned, idents=c("SNU719","NCC24"))
DefaultAssay(coculture.Aza.cells) <- "RNA"
coculture.Aza.cells <- NormalizeData(coculture.Aza.cells)
coculture.Aza.cells <- ScaleData(coculture.Aza.cells)
DefaultAssay(coculture.Aza.cells) <- "SCT"

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
  scale_fill_distiller(palette = "RdBu") + ggtitle("Aza-treated SNU719 and DMSO-treated NCC24 cells in monoculture")

# Plot the same DEG on the mixed.Aza.cells
DoHeatmap(coculture.Aza.cells, features = CellLine.markers.FindAll.heatmap$gene, group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
  scale_fill_distiller(palette = "RdBu") + ggtitle("Aza-treated SNU719 and Aza-treated NCC24 cells in co-culture")
