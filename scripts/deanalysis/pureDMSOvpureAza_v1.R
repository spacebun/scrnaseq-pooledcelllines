#!/usr/bin/env Rscript

# README ------------------------------------------------------------------

# Script version: 0
# Question: How does Aza affect the expression of SNU719 cells for monoculture?

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
setwd("/Users/serena/scrnaseq/results/analysis/pureDMSOvpureAza")

## Load the datasets: All 4 SNU719 datasets, to be integrated into one object
load("/Users/serena/scrnaseq/results/pure_indiv/april2017/SNU719-Aza-april2017/v1/SCT-SNU719-Aza-april2017-v1.RData")
load("/Users/serena/scrnaseq/results/pure_indiv/april2017/SNU719-DMSO-april2017/v1/SCT-SNU719-DMSO-april2017-v1.RData")
load("/Users/serena/scrnaseq/results/pure_indiv/july2017/SNU719-DMSO-july2017/v1/SCT-SNU719-DMSO-july2017-v1.RData")
load("/Users/serena/scrnaseq/results/pure_indiv/sept2017/SNU719-Aza-sept2017/v1/SCT-SNU719-Aza-sept2017-v1.RData")

# Integrate SNU719 Datasets to create one Seurat Object -----------------------------------------------------------
SNU.list <- list(SNU719.DMSO.april2017.seuratobject,SNU719.DMSO.july2017.seuratobject,SNU719.Aza.april2017.seuratobject,SNU719.Aza.sept2017.seuratobject)
SNU.features <- SelectIntegrationFeatures(object.list = SNU.list, nfeatures = 3000)
SNU.list <- PrepSCTIntegration(object.list = SNU.list, anchor.features = SNU.features)
SNU.anchors <- FindIntegrationAnchors(object.list = SNU.list, normalization.method = "SCT", anchor.features = SNU.features)
SNU.integrated <- IntegrateData(anchorset = SNU.anchors, normalization.method = "SCT")

DefaultAssay(SNU.integrated) <- "integrated"
# Save integrated seurat object
save(SNU.integrated, file="/Users/serena/scrnaseq/results/analysis/pureDMSOvpureAza/integrate_SNU.RData")

# Investigate if there are clusters within dataset-----------------------------------------------------------
# Thoughts: If no distinct clusters it may be better to merge with NCC then compare the
# difference in gene expression for DMSO v NCC and Aza v NCC.
# But else tbh should change idents to CellLine then Find All Markers.
#Assess dimensionality
SNU.integrated <- RunPCA(SNU.integrated, verbose = FALSE)
ElbowPlot(SNU.integrated)
# Dims around 15
PCAPlot(SNU.integrated,
        group.by = "sample")
# Run UMAP, FindNeighbors, FindClusters on integrated object for 15 reductions
SNU.integrated <- RunUMAP(SNU.integrated, reduction = "pca", dims = 1:15)

SNU.integrated <- FindNeighbors(SNU.integrated,
                                 dims = 1:15)

SNU.integrated <- FindClusters(SNU.integrated,
                                resolution = c(0.3,0.6,1.0,1.6,2.0),
                                dims = 1:15)

### Explore different resolutions specified in FindClusters
# Assign identity of clusters
Idents(object = SNU.integrated) <- "integrated_snn_res.0.3"
DimPlot(SNU.integrated, reduction = "umap")

# Save integrated seurat object
save(SNU.integrated, file="/Users/serena/scrnaseq/results/analysis/pureDMSOvpureAza/integrate_SNU_withUMAP.RData")

#### ANALYSIS
#### DE Between Aza and DMSO ----------------------------------------------------------

## Set-Up object
# Save back-up
SNU.integrated.backup  <- SNU.integrated 
# Load annotations
annotations <- read.csv("/users/serena/scrnaseq/data/annotation.csv")
# Set assay to RNA
DefaultAssay(SNU.integrated) <- "RNA" 
# Normalize for FindMarkers
SNU.integrated <- NormalizeData(SNU.integrated, verbose = FALSE) 
# FVF to plot 
SNU.integrated <- FindVariableFeatures(SNU.integrated, selection.method = "vst", nfeatures = "3000")
top10 <- head(VariableFeatures(SNU.integrated), 10)
p3 <- VariableFeaturePlot(SNU.integrated)
LabelPoints(plot = p3, points = top10, repel = TRUE)
# ScaleData for Heatmap
SNU.integrated <- ScaleData(SNU.integrated, vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))

# save object
save(SNU.integrated, file = "integrate_SNU_scaled.Rdata")

## FindMarkers between Aza and DMSO
# Set ident to Treatment
Idents(SNU.integrated) <- SNU.integrated$Treatment
# FindMarkers for SNU719 between Aza and DMSO
pure.SNU719.Aza.markers <- FindMarkers(SNU.integrated, ident.1 = "Aza", ident.2 = "DMSO")
# Combine markers with gene descriptions
pure.SNU719.Aza.markers <- pure.SNU719.Aza.markers %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
# List top 20 up reg and top 20 down reg
plot.genes <- (pure.SNU719.Aza.markers %>% top_n(n = 20, wt = avg_logFC))$gene
plot.genes <- append(plot.genes, (pure.SNU719.Aza.markers %>% top_n(n = 20, wt = -avg_logFC))$gene)
# Plot Heatmap of these genes
DoHeatmap(SNU.integrated, features = plot.genes, group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu")

# Average Expression-----------------------------------------------------------
# Plot average expression, and use Hoverlocator to find genes that are outliers
avg.SNU.integrated <- log1p(AverageExpression(SNU.integrated, verbose=FALSE)$RNA)
p1 <- ggplot(avg.SNU.integrated, aes(Aza, DMSO)) + geom_point() + ggtitle("SNU719-DMSO vs SNU719-Aza")
HoverLocator(plot = p1)

# Based on plot, list some outlier genes and display plot with genes labeled
# Can repeat this step with DE genes
outliers <- c("LYZ", "PHGR1", "ID3","IGFBP2","NTM",
              "TFF1","ADIRF", "S100P","AREG")
# outliers <- c("SPINK1", "CST1", "LCN2","ODAM", "PRSS2","REG1A","S100P","TACC1","BPIFB1","PLA2G2A","SLC39A4",
#               "CAMK2N1","CTAG2", "MTA1", "PHGR1", "ALDH3A1", "REG4", "CCK")

p1 <- ggplot(avg.SNU.integrated , aes(DMSO, Aza)) + geom_point() + ggtitle("SNU719-DMSO vs SNU719-Aza")
LabelPoints(p1, points = outliers)

# Feature Plot: Explore the marker genes  -----------------------------------------------------------
# Marker genes: High avglog_FC, pct1 >> pct2.
FeaturePlot(SNU.integrated, features = c("ALDH3A1", "PARK7", "LYZ", "S100P"))

Idents(SNU.integrated) <- "CellLine"
markers.to.plot <- c("ALDH3A1","AGR2", "PARK7", "S100P", "SPINK1","ODAM")
DotPlot(SNU.integrated, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "Treatment") + RotatedAxis()

Idents(SNU.integrated) <- "Treatment"

### Find genes upreg after Aza, for RNA and Integrated Assay----------------------------------------------------------

## Set Assay to RNA
DefaultAssay(SNU.integrated) <- "RNA"
# Find genes up-regulated after Aza
upreg <- FindMarkers(SNU.integrated, ident.1 = "Aza", ident.2 = "DMSO", only.pos = TRUE)
upreg <- upreg %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# List top 20 upreg genes by highest avg_logFC
top20 <- upreg %>% top_n(20, wt = avg_logFC)

# Plot heatmap of these genes
DoHeatmap(SNU.integrated, features = top20$gene, group.bar = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu")
# ----------------------------------------------------------
## Set Assay to integrated
DefaultAssay(SNU.integrated) <- "integrated"
# Find genes up-regulated after Aza
upreg.int <- FindMarkers(SNU.integrated, ident.1 = "Aza", ident.2 = "DMSO", only.pos = TRUE)
upreg.int <- upreg.int %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# List top 100  upreg genes by highest avg_logFC
top100 <- upreg.int %>% top_n(100, wt = avg_logFC)

# Plot heatmap of these genes          
DoHeatmap(SNU.integrated, features = top100$gene, group.bar = TRUE, size = 4, hjust = 1, angle = 0)

# Misc Code -----------------------------------------------------------

# # Plot avg expression of expected DE genes (based on paper)
# aza.genes <- c("p16", "DAPK", "MGMT", "FHIT", "CDKN2B", "ESR1", "IGSF4")
# p2 <- ggplot(avg.SNU.integrated , aes(DMSO, Aza)) + geom_point() + ggtitle("SNU719-DMSO vs SNU719-Aza")
# LabelPoints(p2, points = aza.genes)

# Are these genes in papers?
pure.SNU719.Aza.markers$gene[pure.SNU719.Aza.markers$gene %in% silencedgenes$Gene]
DoHeatmap(SNU.integrated, features = c("TFF1", "CTNNB1"))

# Visualizing DE Genes  -----------------------------------------------------------
# Genes with strong differences based on trial and error
t.de.genes <- c("")

# list top 10 avg log fc
top10 <- top_n(pure.Treatment.markers$gene, 10)
top3 <- top_n(pure.Treatment.markers, 3, avg_logFC)$gene
top20 <- top_n(pure.Treatment.markers, 20, avg_logFC)$gene
# Feature Plot -> test genes based on list^
FeaturePlot(SNU.integrated, 
            features = c(top3), 
            split.by = "Treatment", 
            cols = c("white", "red"),
            label.size = 2)

FeaturePlot(SNU.integrated, 
            features = c("SLCO4A1-AS1"), 
            cols = c("grey", "red"))

# HeatMap 
DoHeatmap(SNU.integrated, features = top20, group.bar = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu")

