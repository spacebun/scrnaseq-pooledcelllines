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

# RunPCA and RunUMAP on integrated object-----------------------------------------------------------
SNU.integrated <- RunPCA(SNU.integrated, verbose = FALSE)
#Assess dimensionality by observing Elbow Plot and PCA Plot
ElbowPlot(SNU.integrated)
PCAPlot(SNU.integrated,
        group.by = "sample")
# Run UMAP, FindNeighbors, FindClusters on integrated object for 15 reductions
SNU.integrated <- RunUMAP(SNU.integrated, reduction = "pca", dims = 1:15)
SNU.integrated <- FindNeighbors(SNU.integrated,dims = 1:15)
SNU.integrated <- FindClusters(SNU.integrated,resolution = c(0.3,0.6,1.0),dims = 1:15)
# Observe UMAP 
Idents(object = SNU.integrated) <- "integrated_snn_res.0.3"
DimPlot(SNU.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# sample
png('Groupby-Sample.png')
DimPlot(SNU.integrated, reduction = "umap", group.by = "sample", label = "TRUE")
dev.off()
# Treatment: Aza / DMSO
png('Groupby-Treatment.png')
DimPlot(SNU.integrated, reduction = "umap", group.by = "Treatment", label = "TRUE")
dev.off()

# Save integrated seurat object
save(SNU.integrated, file="/Users/serena/scrnaseq/results/analysis/pureDMSOvpureAza/integrate_SNU_withUMAP.RData")

# Normalize data in RNA assay-----------------------------------------------------------
DefaultAssay(SNU.integrated) <- "RNA"
SNU.integrated <- NormalizeData(SNU.integrated, verbose = FALSE)
SNU.integrated <- FindVariableFeatures(SNU.integrated, selection.method = "vst", nfeatures = "3000")
top10 <- head(VariableFeatures(SNU.integrated), 10)
p3 <- VariableFeaturePlot(SNU.integrated)
LabelPoints(plot = p3, points = top10, repel = TRUE)
SNU.integrated <- ScaleData(SNU.integrated, vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))

# Set Idents-----------------------------------------------------------
# CellLine.Treatment
Idents(SNU.integrated)<-SNU.integrated$CellLine
SNU.integrated$CellLine.Treatment<- paste(Idents(SNU.integrated), SNU.integrated$Treatment, sep = "_")
# CellLine.Culture
Idents(SNU.integrated)<-SNU.integrated$CellLine
SNU.integrated$CellLine.Culture <- paste(Idents(SNU.integrated), SNU.integrated$Culture, sep = "_")
# CellLine.Treatment.Culture
Idents(SNU.integrated)<-SNU.integrated$CellLine.Treatment
SNU.integrated$CellLine.Treatment.Culture <- paste(Idents(SNU.integrated), SNU.integrated$Culture, sep = "_")
# Save normalized and scaled object with corrected idents
save(SNU.integrated, file="/Users/serena/scrnaseq/results/analysis/pureDMSOvpureAza/integrate_SNU_scaledRNAassay.RData")

# Average Expression-----------------------------------------------------------
DefaultAssay(SNU.integrated) <- "RNA"
Idents(SNU.integrated) <- SNU.integrated$Treatment

# Plot average expression, and use Hoverlocator to find genes that are outliers
avg.SNU.integrated <- log1p(AverageExpression(SNU.integrated, verbose=FALSE)$SCT)
p1 <- ggplot(avg.SNU.integrated, aes(Aza, DMSO)) + geom_point() + ggtitle("SNU719-DMSO vs SNU719-Aza")
HoverLocator(plot = p1)

# DE Analysis-----------------------------------------------------------
DefaultAssay(SNU.integrated) <- "RNA"
Idents(SNU.integrated)<-SNU.integrated$Treatment
# Load annotations
annotations <- read.csv("/users/serena/scrnaseq/data/annotation.csv")

# Find Marker Genes between two Cell Lines
pure.SNU719.Treatment.markers <- FindAllMarkers(SNU.integrated, idents= 'Treatment', only.pos = TRUE)
pure.SNU719.Treatment.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# Combine markers with gene descriptions
pure.SNU719.Treatment.markers <- pure.SNU719.Treatment.markers %>%
  rownames_to_column(var="genes") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# Save list of DE genes
write.csv(pure.SNU719.Treatment.markers,file = "/Users/serena/scrnaseq/results/analysis/pureDMSOvpureAza/All_DE_genes_SNU719_pureDMSOvpureAza.csv")
save(pure.SNU719.Treatment.markers, file = "/Users/serena/scrnaseq/results/analysis/pureDMSOvpureAza/All_DE_genes_SNU719_pureDMSOvpureAza.Rdata")

# Visualizing DE Genes  -----------------------------------------------------------
DefaultAssay(SNU.integrated) <- "RNA"
Idents(SNU.integrated)<-SNU.integrated$Treatment
# subset top 20 DE genes for each treatment condition
top20.pure.SNU719.Treatment.markers <- pure.SNU719.Treatment.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
# plot heatmap of top 20 for each treatment
DoHeatmap(SNU.integrated, features = top20.pure.SNU719.Treatment.markers$gene, 
          group.bar = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu") +
  NoLegend()

# Feature Plot of selected genes across treatment condition
# modify genes under 'features'
FeaturePlot(SNU.integrated, 
            features = c("TFF1"), 
            split.by = "Treatment", 
            cols = c("white", "red"),
            label.size = 2)
# Violin plot of selected genes across treatment condition
# modify genes under 'features'
plots <- VlnPlot(SNU.integrated, 
                 features = c("TMSB4X","AREG","S100P","S100A10","NEAT1", "DCBLD2"), 
                 split.by = "Treatment", 
                 group.by = 'CellLine',
                 pt.size = 0, 
                 combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x +
                  theme(axis.title = element_text(size=8),
                        axis.text.x = element_text(size=8,angle = 0, hjust = 0.5),
                        axis.text.y = element_text(size=8),
                        plot.title=element_text(size=8,face="bold"),
                        legend.position = 'none'))
CombinePlots(plots = plots, ncol = 2)


# OLD CODE ----------------------------------------------------------
# 22 March Identify conserved cell line markers  ----------------------------------------------------------
DefaultAssay(SNU.integrated) <- "RNA"
SNU.integrated <- NormalizeData(SNU.integrated)
SNU.integrated <- FindVariableFeatures(SNU.integrated, selection.method = "vst")
SNU.integrated <- ScaleData(SNU.integrated, features = VariableFeatures(SNU.integrated))

Idents(SNU.integrated) <- "CellLine"
# Identify genes that are conserved markers irrespective of treatment with Aza in the SNU719 and NCC24 clusters
# Cant run?
pure.SNU719.markers <- FindConservedMarkers(SNU.integrated, ident.1 = 'SNU719', grouping.var = "Treatment")
head(pure.SNU719.markersSNU719.markers)

# Explore the marker genes  in FP
# Marker genes: High avglog_FC, pct1 >> pct2.
FeaturePlot(SNU.integrated, features = c("ALDH3A1", "PARK7", "LYZ", "S100P"))

Idents(SNU.integrated) <- "CellLine"
markers.to.plot <- c("ALDH3A1","AGR2", "PARK7", "S100P", "SPINK1","ODAM")
DotPlot(SNU.integrated, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "Treatment") + RotatedAxis()
#  DE Genes between clusters -----------------------------------------------------------
Idents(object = SNU.integrated) <- "integrated_snn_res.0.3"
cluster.markers <- FindAllMarkers(SNU.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(SNU.integrated, features = top10$gene, group.bar = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu")
# Find DE genes with min difff pct -----------------------------------------------------------
DefaultAssay(SNU.integrated) <- "RNA"
pure.Treatment.markers.min <- FindAllMarkers(SNU.integrated, idents= 'Treatment', min.diff.pct = 0.5)
pure.Treatment.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# Combine markers with gene descriptions
pure.Treatment.markers.min_ann <- pure.Treatment.markers.min %>%
  rownames_to_column(var="genes") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
# Displaying DE genes 3 April ----------------------------
Idents(SNU.integrated) <- "Treatment"
top10 <- pure.Treatment.markers_ann %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(SNU.integrated, features = top10$gene, group.bar = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu")
+ NoLegend()
