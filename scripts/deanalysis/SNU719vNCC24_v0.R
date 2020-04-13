#!/usr/bin/env Rscript

# README ------------------------------------------------------------------

# Script version: 0
# Question: How do SNU719 and NCC24 differ in gene expression?
# Using pure SNU-DMSO and pure NCC24-DMSO


# O/P:


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
library(RColorBrewer)

options(future.globals.maxSize = 4000 * 1024^2)

## Set working directory
setwd("/Users/serena/scrnaseq/results/analysis/SNU719vsNCC24")

## Load datasets to be merged.
load("/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/merged_SNUintegrated_NCC_DMSO.RData")

# Transfer object to new object for analysis
snu.ncc <- pure.DMSO.merged

# Examining average expression of genes between datasets-----------------------------------------------------------

# Switch current ident to column holding CellLine metadata
Idents(snu.ncc)<-snu.ncc$CellLine

# Calculate the average gene expression?
avg.snu.ncc <- log1p(AverageExpression(snu.ncc, verbose=FALSE)$SCT)

# Plot average expression, and use Hoverlocator to find genes that are outliers
p1 <- ggplot(avg.snu.ncc, aes(SNU719, NCC24)) + geom_point() + ggtitle("SNU719 vs NCC24")
HoverLocator(plot = p1) 

# Based on plot, list some outlier genes and display plot with genes labeled
# Can repeat this step with DE genes
outliers <- c("SPINK1", "CST1", "LCN2","ODAM", "PRSS2","REG1A","S100P","TACC1","BPIFB1","PLA2G2A","SLC39A4",
              "CAMK2N1","CTAG2", "MTA1", "PHGR1", "ALDH3A1", "REG4", "CCK")
p1 <- ggplot(avg.snu.ncc, aes(SNU719, NCC24)) + geom_point() + ggtitle("SNU719 vs NCC24") 
LabelPoints(p1, points = outliers)

# Use FindMarkers to find DE Genes -----------------------------------------------------------
# Switch current ident to the column holding CellLine Labels.
Idents(snu.ncc)<-snu.ncc$CellLine
annotations <- read.csv("/users/serena/scrnaseq/data/annotation.csv")

# Find Marker Genes between two Cell Lines
CellLine.markers <- FindAllMarkers(snu.ncc, idents= 'CellLine', min.pct = 0.25, logfc.threshold = 0.25)
CellLine.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# Combine markers with gene descriptions
CellLine.markers_ann <- CellLine.markers %>%
  rownames_to_column(var="genes") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# Save before continuing any further
write.csv(CellLine.markers_ann,"Markers for SNU719 vs NCC24.csv")
save(CellLine.markers, file = "/Users/serena/scrnaseq/results/analysis/SNU719vsNCC24/MarkersSNUvNCC.Rdata")
save(CellLine.markers_ann, file = "/Users/serena/scrnaseq/results/analysis/SNU719vsNCC24/MarkersSNUvNCC_ann.Rdata")

# Find Markers again, but adjust min.diff.pct
CellLine.markers.min <- FindAllMarkers(snu.ncc, idents= 'CellLine', min.pct = 0.25, logfc.threshold = 0.25, min.diff.pct = 0.8, )
CellLine.markers.min %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# Combine markers with gene descriptions
CellLine.markers.min_ann <- CellLine.markers.min %>%
  rownames_to_column(var="genes") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

DoHeatmap(snu.ncc, features = head(CellLine.markers.min_ann$gene, n = 10), group.bar = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu")

# Visualizing DE Genes  -----------------------------------------------------------
# Genes with strong differences based on trial and error
de.genes <- c("PHGR1", "ALDH3A1","CAMK2N1", "CTAG2", "GAL", "ATOX1",
  "SPINK1", "ODAM", "S100P","TACC1", "PARK7", "BPIFB1")

# Feature Plot
FeaturePlot(snu.ncc, 
            features = c("GAL", "S100A14", "REG4", "NTM","ATOX1",
                         "S100P","TACC1","BPIFB1","REG1A", "PARK7"), 
            split.by = "CellLine", 
            cols = c("grey", "red"))

FeaturePlot(snu.ncc, 
            features = c("GAL", "S100A14", "REG4", "NTM","ATOX1",
                         "S100P","TACC1","BPIFB1","REG1A", "PARK7"), 
            cols = c("grey", "red"))

# HeatMap for SNU-DMSO merged with NCC-DMSO (pure.DMSO.merged / snu.ncc) 
DoHeatmap(snu.ncc, features = de.genes, group.bar = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu")

# HeatMap for SNU-DMSO merged with NCC-DMSO (pure.DMSO.merged / snu.ncc) 
DoHeatmap(snu.ncc, features = de.genes, group.bar = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu")


# Dot Plot


# Testing the list of de.genes with heatmaps on several other datasets ------------------------------------------------------------------

# HM for SNU-Aza merged with NCC-DMSO (pure.Aza.merged)
load("/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/Aza reference/merged_SNUintegrated_NCC_Aza.RData")
Idents(pure.Aza.merged)<-pure.Aza.merged$CellLine
DoHeatmap(pure.Aza.merged, features = de.genes, group.bar = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu")
  
# HM for SNU-DMSO-Aza-integrated merged with NCC-DMSO (SNU.integrated merged with NCC)
# First load SNU-Integrated object, to be merged with NCC-DMSO
load("/Users/serena/scrnaseq/results/analysis/pureDMSOvpureAza/integrate_SNU_withUMAP.RData")

# Merge SNU-DMSO-Aza with NCC, run standard steps
All.SNU.NCC.merged <- merge(x = NCC24.DMSO.sept2017.seuratobject,
                         y = c(SNU.integrated),
                         add.cell.id = c("NCC24.DMSO.sept2017","SNU719.integrated"),
                         project = "merge_all_SNU_with_NCC",
                         merge.data = TRUE)
All.SNU.NCC.merged <-SCTransform(All.SNU.NCC.merged,
                              return.only.var.genes = F,
                              vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
All.SNU.NCC.merged <- RunPCA(All.SNU.NCC.merged, verbose = FALSE)
All.SNU.NCC.merged <- RunUMAP(All.SNU.NCC.merged, reduction = "pca",dims = 1:15)
All.SNU.NCC.merged <- FindNeighbors(All.SNU.NCC.merged, dims = 1:15)
All.SNU.NCC.merged <- FindClusters(All.SNU.NCC.merged,
                                resolution = c(0.05,0.07,0.1,0.2,0.5))
Idents(object = All.SNU.NCC.merged) <- "SCT_snn_res.0.1"
DimPlot(All.SNU.NCC.merged, reduction = "umap", group.by = "sample")
DimPlot(All.SNU.NCC.merged, reduction = "umap", group.by = "CellLine",
        cols = c('SNU719' = 'tomato2','NCC24' = 'skyblue2'))
DimPlot(All.SNU.NCC.merged, reduction = "umap", group.by = "Treatment",
        cols = c('DMSO' = 'seagreen3','Aza' = 'lightpink1'))
save(All.SNU.NCC.merged, file = "/Users/serena/scrnaseq/results/analysis/SNU719vsNCC24/AllSNUNCCmerged.Rdata")

Idents(All.SNU.NCC.merged)<-All.SNU.NCC.merged$CellLine
DoHeatmap(All.SNU.NCC.merged, features = de.genes, group.bar = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu")

# HM for mixed_DMSO int with mixed_Aza
# Not done yet