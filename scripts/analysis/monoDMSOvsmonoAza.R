#!/usr/bin/env Rscript

# README ------------------------------------------------------------------
# Question: How does Aza affect the expression of SNU719 cells in monoculture?
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
library(RColorBrewer)

options(future.globals.maxSize = 4000 * 1024^2)

## Set working directory
setwd("/scrnaseq/results/analysis/monoDMSOvsmonoAza")

## Load the datasets: All 4 SNU719 datasets, to be integrated into one object
load("/scrnaseq/results/dataset preparation/mono/april2017/SNU719-DMSO-april2017/SCT-SNU719-DMSO-april2017.RData")
load("/scrnaseq/results/dataset preparation/mono/july2017/SNU719-DMSO-july2017/SCT-SNU719-DMSO-july2017.RData")
load("/scrnaseq/results/dataset preparation/mono/april2017/SNU719-Aza-april2017/SCT-SNU719-Aza-april2017.RData")
load("/scrnaseq/results/dataset preparation/mono/sept2017/SNU719-Aza-sept2017/SCT-SNU719-Aza-sept2017.RData")

# Integrate SNU719 Datasets to create one Seurat Object -----------------------------------------------------------
SNU.list <- list(SNU719.DMSO.april2017.seuratobject,SNU719.DMSO.july2017.seuratobject,SNU719.Aza.april2017.seuratobject,SNU719.Aza.sept2017.seuratobject)
SNU.features <- SelectIntegrationFeatures(object.list = SNU.list, nfeatures = 3000)
SNU.list <- PrepSCTIntegration(object.list = SNU.list, anchor.features = SNU.features)
SNU.anchors <- FindIntegrationAnchors(object.list = SNU.list, normalization.method = "SCT", anchor.features = SNU.features)
SNU.integrated <- IntegrateData(anchorset = SNU.anchors, normalization.method = "SCT")
DefaultAssay(SNU.integrated) <- "integrated"

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
# Save integrated seurat object, normalized and scaled with corrected idents
save(SNU.integrated, file="integrate_SNU.RData")

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
annotations <- read.csv("/scrnaseq/data/annotation.csv")

# Find Marker Genes between two Cell Lines
pure.SNU719.Treatment.markers <- FindAllMarkers(SNU.integrated, idents= 'Treatment', only.pos = TRUE)
pure.SNU719.Treatment.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# Combine markers with gene descriptions
pure.SNU719.Treatment.markers <- pure.SNU719.Treatment.markers %>%
  rownames_to_column(var="genes") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

length(which(pure.SNU719.Treatment.markers$cluster=="SNU719_Aza"))

# Save list of DE genes
write.csv(pure.SNU719.Treatment.markers,file = "DEG_monoDvsmonoA.csv")
save(pure.SNU719.Treatment.markers, file = "DEG_monoDvsmonoA.Rdata")

# Visualizing DE Genes  -----------------------------------------------------------
DefaultAssay(SNU.integrated) <- "RNA"
Idents(SNU.integrated)<-SNU.integrated$Treatment
# subset top 20 DE genes for each treatment condition
top20.pure.SNU719.Treatment.markers <- pure.SNU719.Treatment.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
# plot heatmap of top 20 for each treatment

DoHeatmap(SNU.integrated, features = top20.pure.SNU719.Treatment.markers$gene,
          group.bar = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu")

# replot absolute values
DoHeatmap(SNU.integrated, features = top20.pure.SNU719.Treatment.markers$gene,
          group.bar = TRUE, size = 4, hjust = 1, angle = 0,
          slot = "counts") +
  scale_fill_distiller(palette = "RdBu")

# Violin plot of selected genes across treatment condition
# modify genes under 'features'
plots <- VlnPlot(SNU.integrated,
                 features = c("TMSB4X","AREG","S100P","S100A10","NEAT1", "DCBLD2"),
                 split.by = "Treatment",
                 group.by = 'CellLine',
                 pt.size = 0,
                 combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x +
                  theme(axis.title = element_text(size=15),
                        axis.text.x = element_text(size=15,angle = 0, hjust = 0.5),
                        axis.text.y = element_text(size=15),
                        plot.title=element_text(size=20,face="bold"),
                        legend.position = 'top'))
CombinePlots(plots = plots, ncol = 2)


# FeaturePlot(SNU.integrated, features = c("DCBLD2"), split.by = "CellLine.Treatment", sort.cell = TRUE)
# DotPlot(SNU.integrated, features = c("DCBLD2")) + RotatedAxis()
# DoHeatmap(subset(SNU.integrated, downsample = 100), features =  c("TMSB4X","AREG","S100P","S100A10","NEAT1", "DCBLD2"), size = 3)
