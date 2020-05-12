#!/usr/bin/env Rscript

# README ------------------------------------------------------------------
# Question: How do culture conditions affect the expression of DMSO-treated SNU719 cells?

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
setwd("/scrnaseq/results/analysis/monoDMSOvscoDMSO")

## Load the datasets: 2 SNU719-DMSO, 1 mixed_DMSO_labeled
load("/scrnaseq/results/dataset preparation/mono/april2017/SNU719-DMSO-april2017/SCT-SNU719-DMSO-april2017.RData")
load("/scrnaseq/results/dataset preparation/mono/july2017/SNU719-DMSO-july2017/SCT-SNU719-DMSO-july2017.RData")
load("/scrnaseq/results/cell line assignment/DMSO/coculture-DMSO-assigned.RData")

# Integrate datasets: SNU719 and coculture_DMSO--------------------------------------
DMSO.list <- list(SNU719.DMSO.april2017.seuratobject,SNU719.DMSO.july2017.seuratobject, coculture.DMSO.assigned)
DMSO.features <- SelectIntegrationFeatures(object.list = DMSO.list, nfeatures = 3000)
DMSO.list <- PrepSCTIntegration(object.list = DMSO.list, anchor.features = DMSO.features)
DMSO.anchors <- FindIntegrationAnchors(object.list = DMSO.list, normalization.method = "SCT", anchor.features = DMSO.features)
DMSO.integrated <- IntegrateData(anchorset = DMSO.anchors, normalization.method = "SCT")
DefaultAssay(DMSO.integrated) <- "integrated"

# RunPCA Elbow Plot UMAP
DMSO.integrated <- RunPCA(DMSO.integrated, verbose = FALSE)
ElbowPlot(DMSO.integrated)
PCAPlot(DMSO.integrated,group.by = "sample")
DMSO.integrated <- RunUMAP(DMSO.integrated, reduction = "pca", dims = 1:15)

DimPlot(DMSO.integrated, reduction = "umap") + ggtitle("Integrated Mono and Co SNU719 DMSO Datasets")

DimPlot(DMSO.integrated, group.by = "CellLine.Culture",
        cols = c('SNU719_Co' = 'tomato2', 'SNU719_Mono' = 'lightpink1',
                 'NCC24_Co' = "gray70",'1_Co' = "gray70",'2_Co' = "gray70",'4_Co' = "gray70",'5_Co' = "gray70",'6_Co' = "gray70",'7_Co' = "gray70")) +
  ggtitle("Integrated Mono and Co SNU719 DMSO Datasets")

DimPlot(DMSO.integrated, group.by = "CellLine",
        cols = c('SNU719' = 'tomato2',
                 'NCC24' = "forestgreen",'1' = "gray70",'2' = "gray70",'4' = "gray70",'5' = "gray70",'6' = "gray70",'7' = "gray70")) +
  ggtitle("Integrated Mono and Co SNU719 DMSO Datasets")

Idents(DMSO.integrated) <- "CellLine"
DimPlot(DMSO.integrated, split.by = "Culture")


# Set-up SNU object for DE: Normalize, Set Idents, ScaleData  ----------------------------------------------------------
# Set assay to RNA
DefaultAssay(DMSO.integrated) <- "RNA"
# Normalize RNA data for visualization purposes
DMSO.integrated <- NormalizeData(DMSO.integrated, verbose = FALSE)
# Scale Data, and save the object again
# FVF to plot
DMSO.integrated <- FindVariableFeatures(DMSO.integrated, selection.method = "vst", nfeatures = "3000")
top10 <- head(VariableFeatures(DMSO.integrated), 10)
p3 <- VariableFeaturePlot(DMSO.integrated)
LabelPoints(plot = p3, points = top10, repel = TRUE)
DMSO.integrated <- ScaleData(DMSO.integrated,
                                       features = VariableFeatures(DMSO.integrated),
                                       vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
# Save normalized, scaled object
save(DMSO.integrated, file="DMSO_integrated.RData")

# FCM: SNU Across Culture--------------------------------------
DefaultAssay(SNU719.Culture.Cells) <- "RNA"
Idents(SNU719.Culture.Cells) <- "CellLine"
# #Find conserved genes for SNU719-DMSO across culture (Co v Mono)
SNU719.markers <- FindConservedMarkers(SNU719.Culture.Cells, ident.1 = 'SNU719', grouping.var = "Culture", verbose = FALSE)
head(SNU719.markers)

# Explore the marker genes for  SNU719 in FP
# Marker genes: High avglog_FC, pct1 >> pct2.
FeaturePlot(DMSO.integrated, features = c("SPINK4"))

Idents(DMSO.integrated) <- "CellLine"
markers.to.plot <- c("HIST1H4C","CCK","MT2A","LYZ")
DotPlot(DMSO.integrated, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8,
        split.by = "Culture") + RotatedAxis()


# Average Expression of SNU719-----------------------------------------------------------
Idents(DMSO.integrated) <- DMSO.integrated$CellLine
SNU719.cells <- subset(DMSO.integrated, idents = c("SNU719"))
Idents(SNU719.cells) <-DMSO.integrated$Culture
avg.SNU719.cells <- log1p(AverageExpression(SNU719.cells, verbose = FALSE)$RNA)
p1 <- ggplot(avg.SNU719.cells, aes(`Mono`, `Co`)) + geom_point() + ggtitle("DMSO-treated SNU719 in Mono vs DMSO-treated SNU719 in Co-Culture")
HoverLocator(plot = `p1`)
genes.to.label = c("IFI6", "IFI27", "STAT1", "NOP53", "SLC25A6", "PRDX11","ACB","TPD52L1")
LabelPoints(p1, points = genes.to.label)
#DE: SNU--------------------------------------
# Load annotations
annotations <- read.csv("/scrnaseq/data/annotation.csv")
# Create new object to perform analysis on:
DMSO.integrated.new <- DMSO.integrated
Idents(DMSO.integrated.new)<-DMSO.integrated.new$CellLine.Culture
DefaultAssay(DMSO.integrated.new) <- "RNA"

# FindMarkers
DMSO.SNU719.Culture.response <- FindMarkers(DMSO.integrated.new, ident.1 = "SNU719_Co", ident.2 = "SNU719_Mono", verbose = TRUE, min.diff.pct = 0.3)
DMSO.SNU719.Culture.response <- DMSO.SNU719.Culture.response %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

length(which(DMSO.SNU719.Culture.response$avg_logFC>0))


# save list of de genes
write.csv(DMSO.SNU719.Culture.response,file = "DEG_monoDMSOvscoDMSO.csv")
save(DMSO.SNU719.Culture.response, file = "DEG_monoDMSOvscoDMSO.Rdata")

# Visualizing DE Genes  -----------------------------------------------------------
# List top 20 up reg and top 20 down reg
top20.DMSO.SNU719.Culture.markers <- (top_n(DMSO.SNU719.Culture.response, n = 20, wt = avg_logFC))$gene
top20.DMSO.SNU719.Culture.markers <- append(top20.DMSO.SNU719.Culture.markers, (top_n(DMSO.SNU719.Culture.response, n = 20, wt = -avg_logFC))$gene)

# Create object of just SNU719_Mono and SNU719_Culture cells
Idents(DMSO.integrated)<-DMSO.integrated$CellLine.Culture
DefaultAssay(DMSO.integrated) <- "RNA"
DMSO.SNU719.Cells <- subset(DMSO.integrated, idents = c("SNU719_Mono", "SNU719_Co"))

# Plot Heatmap of these genes
DefaultAssay(DMSO.SNU719.Cells) <- "RNA"
DoHeatmap(DMSO.SNU719.Cells, group.by = "CellLine.Culture",features = top20.DMSO.SNU719.Culture.markers,
          group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
  scale_fill_distiller(palette = "RdBu")

# replot absolute values
DoHeatmap(DMSO.SNU719.Cells, features = top20.DMSO.SNU719.Culture.markers,
          group.bar = TRUE, size = 4, hjust = 1, angle = 0,
          slot = "counts") +
  scale_fill_distiller(palette = "RdBu")

# Feature Plot of selected genes across treatment condition
# modify genes under 'features'
# moddify object to be plotted as desiredzw
FeaturePlot(DMSO.SNU719.Cells, features = c("CDKN2A","ID2","IGFBP3"),
            split.by = "Treatment", max.cutoff = 3, cols = c("grey", "red"))

# Violin plot of selected genes across treatment condition
# modify genes under 'features'
plots <- VlnPlot(DMSO.SNU719.Cells,
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
DMSO.SNU719.Cells.test <-DMSO.SNU719.Cells
DMSO.SNU719.Cells.test <-DMSO.SNU719.Cells
DMSO.SNU719.Cells.test.subset <-subset(DMSO.SNU719.Cells.test, subset = IFI27 >1)

plots <- VlnPlot(DMSO.SNU719.Cells.test.subset,
                 features = c("IFI27","CD9","IFI6", "CD81"),
                 split.by = "Culture",
                 group.by = 'CellLine',
                 pt.size = 0,
                 combine = FALSE,
                 # slot = "counts"
)
plots <- lapply(X = plots, FUN = function(x) x +
                  # scale_y_continuous(breaks=seq(0,13,1)) +
                  theme(axis.title = element_text(size=8),
                        axis.text.x = element_text(size=8,angle = 0, hjust = 0.5),
                        axis.text.y = element_text(size=8),
                        plot.title=element_text(size=8,face="bold"),
                        legend.position = 'top'))
CombinePlots(plots = plots, nrow = 2)
