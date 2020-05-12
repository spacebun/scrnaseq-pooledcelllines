#!/usr/bin/env Rscript

# README ------------------------------------------------------------------
# Question: How do culture conditions affect the expression of Aza-treated SNU719 cells?

# Set-up packages, working directory, datasets, set idents -----------------------------------------------------------
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
setwd("/scrnaseq/results/analysis/monoAzavscoAza")

## Load the datasets: 2 SNU719-DMSO, 1 NCC24-DMSO,1 mixed_DMSO_labeled
load("/scrnaseq/results/dataset preparation/mono/april2017/SNU719-Aza-april2017/SCT-SNU719-Aza-april2017.RData")
load("/scrnaseq/results/dataset preparation/mono/sept2017/SNU719-Aza-sept2017/SCT-SNU719-Aza-sept2017.RData")
load("/scrnaseq/results/cell line assignment/Aza/coculture-Aza-assigned.RData")

# Integrate datasets: SNU719 and mixed_Aza--------------------------------------
Aza.list <- list(SNU719.Aza.april2017.seuratobject,SNU719.Aza.sept2017.seuratobject, coculture.Aza.assigned)
Aza.features <- SelectIntegrationFeatures(object.list = Aza.list, nfeatures = 3000)
Aza.list <- PrepSCTIntegration(object.list = Aza.list, anchor.features = Aza.features)
Aza.anchors <- FindIntegrationAnchors(object.list = Aza.list, normalization.method = "SCT", anchor.features = Aza.features)
Aza.integrated <- IntegrateData(anchorset = Aza.anchors, normalization.method = "SCT")
DefaultAssay(Aza.integrated) <- "integrated"

# RunPCA Elbow Plot UMAP
Aza.integrated <- RunPCA(Aza.integrated, verbose = FALSE)
ElbowPlot(Aza.integrated)
PCAPlot(Aza.integrated,group.by = "sample")
Aza.integrated <- RunUMAP(Aza.integrated, reduction = "pca", dims = 1:15)

DimPlot(Aza.integrated, reduction = "umap") + ggtitle("Integrated Mono and Co SNU719 Aza Datasets")

DimPlot(Aza.integrated, group.by = "CellLine.Culture",
        cols = c('SNU719_Co' = 'tomato2', 'SNU719_Mono' = 'lightpink1',
                 'NCC24_Co' = "gray70",'1_Co' = "gray70",'2_Co' = "gray70",'4_Co' = "gray70",'5_Co' = "gray70",'6_Co' = "gray70",'7_Co' = "gray70")) +
  ggtitle("Integrated Mono and Co SNU719 Aza Datasets")

DimPlot(Aza.integrated, group.by = "CellLine",
        cols = c('SNU719' = 'tomato2',
                 'NCC24' = "forestgreen",'1' = "gray70",'2' = "gray70",'4' = "gray70",'5' = "gray70",'6' = "gray70",'7' = "gray70")) +
  ggtitle("Integrated Mono and Co SNU719 Aza Datasets")

Idents(Aza.integrated) <- "CellLine"
DimPlot(Aza.integrated, split.by = "Culture")


# Set-up SNU object for DE: Normalize, Set Idents, ScaleData  ----------------------------------------------------------
# Set assay to RNA
DefaultAssay(Aza.integrated) <- "RNA"
# Normalize RNA data for visualization purposes
Aza.integrated <- NormalizeData(Aza.integrated, verbose = FALSE)
# Scale Data, and save the object again
# FVF to plot
Aza.integrated <- FindVariableFeatures(Aza.integrated, selection.method = "vst", nfeatures = "3000")
top10 <- head(VariableFeatures(Aza.integrated), 10)
p3 <- VariableFeaturePlot(Aza.integrated)
LabelPoints(plot = p3, points = top10, repel = TRUE)
Aza.integrated <- ScaleData(Aza.integrated,
                             features = VariableFeatures(Aza.integrated),
                             vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
# Save normalized, scaled object
save(Aza.integrated, file="Aza_integrated.RData")

# FCM: SNU Across Culture--------------------------------------
DefaultAssay(SNU719.Culture.Cells) <- "RNA"
Idents(SNU719.Culture.Cells) <- "CellLine"
# #Find conserved genes for SNU719-Aza across culture (Co v Mono)
SNU719.markers <- FindConservedMarkers(SNU719.Culture.Cells, ident.1 = 'SNU719', grouping.var = "Culture", verbose = FALSE)
head(SNU719.markers)

# Explore the marker genes for  SNU719 in FP
# Marker genes: High avglog_FC, pct1 >> pct2.
FeaturePlot(Aza.integrated, features = c("SPINK4"))

Idents(Aza.integrated) <- "CellLine"
markers.to.plot <- c("HIST1H4C","CCK","MT2A","LYZ")
DotPlot(Aza.integrated, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8,
        split.by = "Culture") + RotatedAxis()


# Average Expression of SNU719-----------------------------------------------------------
Idents(Aza.integrated) <- Aza.integrated$CellLine
SNU719.cells <- subset(Aza.integrated, idents = c("SNU719"))
Idents(SNU719.cells) <-Aza.integrated$Culture
avg.SNU719.cells <- log1p(AverageExpression(SNU719.cells, verbose = FALSE)$RNA)
p1 <- ggplot(avg.SNU719.cells, aes(`Mono`, `Co`)) + geom_point() + ggtitle("Aza-treated SNU719 in Mono vs Aza-treated SNU719 in Co-Culture")
HoverLocator(plot = `p1`)
genes.to.label = c("IFI6", "IFI27", "STAT1", "NOP53", "SLC25A6", "PRDX11","ACB","TPD52L1")
LabelPoints(p1, points = genes.to.label)
#DE: SNU--------------------------------------
# Load annotations
annotations <- read.csv("/scrnaseq/data/annotation.csv")
# Create new object to perform analysis on:
Aza.integrated.new <- Aza.integrated
Idents(Aza.integrated.new)<-Aza.integrated.new$CellLine.Culture
DefaultAssay(Aza.integrated.new) <- "RNA"

# FindMarkers
Aza.SNU719.Culture.response <- FindMarkers(Aza.integrated.new, ident.1 = "SNU719_Co", ident.2 = "SNU719_Mono", verbose = TRUE, min.diff.pct = 0.3)
Aza.SNU719.Culture.response <- Aza.SNU719.Culture.response %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

length(which(Aza.SNU719.Culture.response$avg_logFC>0))

# save list of de genes
write.csv(Aza.SNU719.Culture.response,file = "DEG_monoAzavcoAza.csv")
save(Aza.SNU719.Culture.response, file = "DEG_monoAzavcoAza.Rdata")

# Visualizing DE Genes  -----------------------------------------------------------
# List top 20 up reg and top 20 down reg
top20.Aza.SNU719.Culture.markers <- (top_n(Aza.SNU719.Culture.response, n = 20, wt = avg_logFC))$gene
top20.Aza.SNU719.Culture.markers <- append(top20.Aza.SNU719.Culture.markers, (top_n(Aza.SNU719.Culture.response, n = 20, wt = -avg_logFC))$gene)

# Create object of just SNU719_Mono and SNU719_Culture cells
Aza.SNU719.Cells <- subset(Aza.integrated.new, idents = c("SNU719_Mono", "SNU719_Co"))

# Plot Heatmap of these genes
DefaultAssay(Aza.SNU719.Cells) <- "RNA"
DoHeatmap(Aza.SNU719.Cells, group.by = "CellLine.Culture",features = top20.Aza.SNU719.Culture.markers,
          group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
  scale_fill_distiller(palette = "RdBu")

# replot absolute values
DoHeatmap(Aza.SNU719.Cells, features = top20.Aza.SNU719.Culture.markers,
          group.bar = TRUE, size = 4, hjust = 1, angle = 0,
          slot = "counts") +
  scale_fill_distiller(palette = "RdBu")

# Feature Plot of selected genes across treatment condition
# modify genes under 'features'
# moddify object to be plotted as desiredzw
FeaturePlot(Aza.SNU719.Cells, features = c("CDKN2A","ID2","IGFBP3"),
            split.by = "Treatment", max.cutoff = 3, cols = c("grey", "red"))

# Violin plot of selected genes across treatment condition
# modify genes under 'features'
plots <- VlnPlot(Aza.SNU719.Cells,
                 features = c("IFI27","CD9","IFI6", "CD81" ),
                 split.by = "Culture",
                 group.by = 'CellLine',
                 pt.size = 0,
                 combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x +
                  theme(axis.title = element_text(size=15),
                        axis.text.x = element_text(size=15,angle = 0, hjust = 0.5),
                        axis.text.y = element_text(size=15),
                        plot.title=element_text(size=20,face="bold"),
                        legend.position = 'top'))
CombinePlots(plots = plots, nrow = 2)

# Investigating the outlier in Vln Plots
Aza.SNU719.Cells.test <-Aza.SNU719.Cells
Aza.SNU719.Cells.test <-Aza.SNU719.Cells
Aza.SNU719.Cells.test.subset <-subset(Aza.SNU719.Cells.test, subset = IFI27 >1)

plots <- VlnPlot(Aza.SNU719.Cells.test.subset,
                 features = c("IFI27","CD9","IFI6", "CD81"),
                 split.by = "Culture",
                 group.by = 'CellLine',
                 pt.size = 0,
                 combine = FALSE,
                 # slot = "counts"
)
plots <- lapply(X = plots, FUN = function(x) x +
                  theme(axis.title = element_text(size=15),
                        axis.text.x = element_text(size=15,angle = 0, hjust = 0.5),
                        axis.text.y = element_text(size=15),
                        plot.title=element_text(size=20,face="bold"),
                        legend.position = 'top'))
CombinePlots(plots = plots, nrow = 2)
