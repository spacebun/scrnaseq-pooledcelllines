#!/usr/bin/env Rscript

# README ------------------------------------------------------------------
# Question: How does Aza affect the expression of SNU719 cells for co-culture?
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

setwd("/scrnaseq/results/analysis/coDMSOvscoAza")

load("/scrnaseq/results/cell line assignment/DMSO/coculture-DMSO-assigned.RData")
load("/scrnaseq/results/cell line assignment/Aza/coculture-Aza-assigned.RData")

# Integrate coculture datasets ------------------------------------------------------------------
Idents(coculture.DMSO.assigned) <- coculture.DMSO.assigned$CellLine
Idents(coculture.Aza.assigned) <- coculture.Aza.assigned$CellLine
coculture.list <- list(coculture.DMSO.assigned,coculture.Aza.assigned)
coculture.features <- SelectIntegrationFeatures(object.list = coculture.list, nfeatures = 3000)
coculture.list <- PrepSCTIntegration(object.list = coculture.list, anchor.features = coculture.features)
coculture.anchors <- FindIntegrationAnchors(object.list = coculture.list, normalization.method = "SCT", anchor.features = coculture.features)
coculture.integrated <- IntegrateData(anchorset = coculture.anchors, normalization.method = "SCT")
DefaultAssay(coculture.integrated) <- "integrated"

# RunPCA and RunUMAP on integrated object-----------------------------------------------------------
coculture.integrated <- RunPCA(coculture.integrated, verbose = FALSE)
#Assess dimensionality by observing Elbow Plot and PCA Plot
ElbowPlot(coculture.integrated)
PCAPlot(coculture.integrated,group.by = "sample")
# Run UMAP, FindNeighbors, FindClusters on integrated object for 15 reductions
coculture.integrated <- RunUMAP(coculture.integrated, reduction = "pca", dims = 1:15)
coculture.integrated <- FindNeighbors(coculture.integrated, dims = 1:15)
coculture.integrated <- FindClusters(coculture.integrated,dims = 1:15, resolution = c(0.2,0.3,0.6))
Idents(object = coculture.integrated) <- "integrated_snn_res.0.3"

DimPlot(coculture.integrated, reduction = "umap") + ggtitle("Integrated Co-Cultured Datasets")
DimPlot(coculture.integrated, group.by = "CellLine", split.by = "Treatment", ncol = 2)

# save(coculture.integrated, file="integrate_coculture_withUMAP.RData")

### Set Idents -----------------------------------------------------------
# CellLine.Treatment
Idents(coculture.integrated)<- coculture.integrated$CellLine
coculture.integrated$CellLine.Treatment<- paste(Idents(coculture.integrated), coculture.integrated$Treatment, sep = "_")
# CellLine.Culture
Idents(coculture.integrated)<-coculture.integrated$CellLine
coculture.integrated$CellLine.Culture <- paste(Idents(coculture.integrated), coculture.integrated$Culture, sep = "_")
# CellLine.Treatment.Culture
Idents(coculture.integrated)<-coculture.integrated$CellLine.Treatment
coculture.integrated$CellLine.Treatment.Culture <- paste(Idents(coculture.integrated), coculture.integrated$Culture, sep = "_")

DimPlot(coculture.integrated, reduction = "umap", group.by = "CellLine.Treatment",
        cols = c('SNU719_Aza' = 'tomato2','SNU719_DMSO' = 'lightpink1','NCC24_Aza' = 'forestgreen', 'NCC24_DMSO' = 'seagreen3',
                 '1_Aza' = 'gray70', '1_DMSO' = 'gray70', '2_Aza' = 'gray70','2_DMSO' = 'gray70', '4_Aza' = 'gray70', '4_DMSO' = 'gray70', '5_Aza' = 'gray70', '5_DMSO' = 'gray70','6_Aza' = 'gray70', '6_DMSO' = 'gray70','7_Aza' = 'gray70', '7_DMSO' = 'gray70')) +
  ggtitle("Integrated Co-Cultured Datasets")
DimPlot(coculture.integrated, group.by = "CellLine", split.by = "Treatment", ncol = 2)
DimPlot(coculture.integrated, group.by = "CellLine.Treatment")
# Save integrated seurat object
save(coculture.integrated, file="integrate_coculture_withUMAP_withIdents.RData")


# Normalize data in RNA assay ----------------------------------------------------------
DefaultAssay(coculture.integrated) <- "RNA"
Idents(coculture.integrated)<-coculture.integrated$CellLine.Treatment
coculture.integrated <- NormalizeData(coculture.integrated, verbose = FALSE)
coculture.integrated <- FindVariableFeatures(coculture.integrated, selection.method = "vst", nfeatures = "3000")
top10 <- head(VariableFeatures(coculture.integrated), 10)
p3 <- VariableFeaturePlot(coculture.integrated)
LabelPoints(plot = p3, points = top10, repel = TRUE)
coculture.integrated <- ScaleData(coculture.integrated, features = VariableFeatures(coculture.integrated), vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
save(coculture.integrated, file="integrate_coculture.RData")


# Average Expression-----------------------------------------------------------
Idents(mixed.integrated) <- mixed.integrated$CellLine
SNU719.cells <- subset(mixed.integrated, idents = c("SNU719"))
Idents(SNU719.cells) <- mixed.integrated$Treatment
avg.SNU719.cells <- log1p(AverageExpression(SNU719.cells, verbose = FALSE)$RNA)
# p1 <- ggplot(avg.SNU719.cells, aes(`DMSO`, `Aza`)) + geom_point() + ggtitle("SNU719-DMSO vs SNU719-Aza in Co-Culture")
# HoverLocator(plot = `p1`)

NCC24.cells <- subset(mixed.integrated, idents = c("NCC24"))
Idents(NCC24.cells) <- "Treatment"
# avg.NCC24.cells <- log1p(AverageExpression(NCC24.cells, verbose = FALSE)$RNA)
# p2 <- ggplot(avg.NCC24.cells, aes(`DMSO`, `Aza`)) + geom_point() + ggtitle("NCC24-DMSO vs NCC24-Aza in Co-Culture")
# HoverLocator(plot = `p2`)
# genes.to.label = c("DKK1","PAGE4","AREG", "CAMK2N1","NFIA", "CCL20", "CXCL10", "TFF1", "CTNNB1","LCN1", "HIST1H4C","SDC4")
# LabelPoints(p2, points = genes.to.label)

#save objects of SNU and NCC
save(SNU719.cells, file = "co-SNU719-cells.Rdata")
save(NCC24.cells, file = "co-NCC24-cells.Rdata")

# DE Analaysis: Comparing DMSO and Aza for SNU719 and NCC24 each ------------------------------------------------------------------
# Load annotations
annotations <- read.csv("/scrnaseq/data/annotation.csv")
# Create new object to perform analysis on:
coculture.integrated.new <- coculture.integrated
Idents(coculture.integrated.new)<-coculture.integrated.new$CellLine.Treatment
DefaultAssay(coculture.integrated.new) <- "RNA"

### SNU719
# FindMarkers
coculture.SNU719.Treatment.markers  <- FindMarkers(coculture.integrated.new, ident.1 = "SNU719_Aza", ident.2 = "SNU719_DMSO", verbose = FALSE)
coculture.SNU719.Treatment.markers <- coculture.SNU719.Treatment.markers %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

length(which(coculture.SNU719.Treatment.markers$avg_logFC<0))

# save list of de genes
write.csv(coculture.SNU719.Treatment.markers,file = "DEG_coDvcoA.csv")
save(coculture.SNU719.Treatment.markers, file = "DEG_coDvcoA.Rdata")
# Visualizing DE Genes  -----------------------------------------------------------
Idents(coculture.integrated.new)<-coculture.integrated.new$CellLine.Treatment
DefaultAssay(coculture.integrated.new) <- "RNA"

### SNU719
# subset top 20 DE genes for each treatment condition
top20.coculture.SNU719.Treatment.markers  <- (top_n(coculture.SNU719.Treatment.markers, n = 20, wt = avg_logFC))$gene
top20.coculture.SNU719.Treatment.markers <- append(top20.coculture.SNU719.Treatment.markers, (top_n(coculture.SNU719.Treatment.markers, n = 20, wt = -avg_logFC))$gene)
# Create object of just SNU719_DMSO and SNU719_Aza cells
Idents(coculture.integrated.new)<-coculture.integrated.new$CellLine.Treatment
coculture.SNU719 <- subset(coculture.integrated.new, idents = c("SNU719_DMSO", "SNU719_Aza"))
Idents(coculture.SNU719)<-coculture.SNU719$CellLine.Treatment
levels(coculture.SNU719)
levels(coculture.SNU719) <- c("SNU719_DMSO","SNU719_Aza")
levels(coculture.SNU719)

# plot heatmap of top 20 for each treatment
DoHeatmap(coculture.SNU719, group.by = "CellLine.Treatment",features = top20.coculture.SNU719.Treatment.markers,
          group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
  scale_fill_distiller(palette = "RdBu")

# replot absolute values
DoHeatmap(coculture.SNU719, features = top20.coculture.SNU719.Treatment.markers,
          group.bar = TRUE, size = 4, hjust = 1, angle = 0,
          slot = "counts") +
  scale_fill_distiller(palette = "RdBu")

# Feature Plot of selected genes across treatment condition
# modify genes under 'features'
# modify object to be plotted as desired
FeaturePlot(coculture.SNU719, features = c("IFI27","IFI6"," LGALS1","RPL17"),
            split.by = "Treatment", max.cutoff = 3, cols = c("grey", "red"))

# Violin plot of selected genes across treatment condition
# modify genes under 'features'
plots <- VlnPlot(coculture.SNU719,
                 features = c("IFI27","IFI6","LGALS1", "RPL17"),
                 split.by = "CellLine.Treatment",
                 group.by = 'CellLine.Treatment',
                 pt.size = 0.5,
                 combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x +
                  theme(axis.title = element_text(size=12),
                        axis.text.x = element_text(size=10,angle = 0, hjust = 0.5),
                        axis.text.y = element_text(size=10),
                        plot.title=element_text(size=8,face="bold"),
                        legend.position = 'top'))
CombinePlots(plots = plots, nrow = 2)

plots <- VlnPlot(coculture.SNU719,
                 features = c("IFI27","IFI6","LGALS1", "RPL17"),
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
CombinePlots(plots = plots, nrow = 2)

# Investigating the bimodality in Vln Plots
coculture.SNU719.test <-coculture.SNU719

coculture.SNU719.test.subset <-subset(coculture.SNU719.test, subset = IFI27 <1)

plots <- VlnPlot(coculture.SNU719.test.subset,
                 features = c("IFI27","IFI6","LGALS1", "RPL17"),
                 split.by = "Treatment",
                 group.by = 'CellLine.Treatment',
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
