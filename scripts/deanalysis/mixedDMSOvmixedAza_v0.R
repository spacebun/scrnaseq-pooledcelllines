#!/usr/bin/env Rscript

# README ------------------------------------------------------------------

# Script version: 0
# Question: How does Aza affect the expression of SNU719 cells for co-culture?
# Compare the two labeled mixed datasets, compare across treatment condition.

# When checking for conserved genes across conditions you can use the integrated assy
# else you should use the RNA assay to compare across conditions.

# O/P:

# Summary of steps
# 1) Integrate the two labeled mixed datasets , run UMAP, label idents -> mixed.integrated
# /Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza/integrate_mixed_withUMAP_updated_Idents.RData
# 2) Which assay to use? Before Testing for differential genes between conditions. (19 March)


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

setwd("/Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza")

# Load mixed objects with cells labeled
load("/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/mixed_DMSO_labeled_v0.RData")
load("/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/mixed_Aza_labeled_v0.RData")

# # Plot the UMAP of both datasets 
# pDMSO <- DimPlot(mixed_DMSO_labeled, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) + 
#   NoLegend() + 
#   ggtitle("Co-Cultured DMSO Dataset")
# pAza <- DimPlot(mixed_Aza_labeled, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) + NoLegend() + ggtitle("Co-Cultured Aza Dataset")
# plot_grid(pDMSO,pAza)

# Integrate mixed datasets ------------------------------------------------------------------
Idents(mixed_DMSO_labeled) <- mixed_DMSO_labeled$CellLine
Idents(mixed_Aza_labeled) <- mixed_Aza_labeled$CellLine
mixed.list <- list(mixed_DMSO_labeled,mixed_Aza_labeled)
mixed.features <- SelectIntegrationFeatures(object.list = mixed.list, nfeatures = 3000)
mixed.list <- PrepSCTIntegration(object.list = mixed.list, anchor.features = mixed.features)
mixed.anchors <- FindIntegrationAnchors(object.list = mixed.list, normalization.method = "SCT", anchor.features = mixed.features)
mixed.integrated <- IntegrateData(anchorset = mixed.anchors, normalization.method = "SCT")
DefaultAssay(mixed.integrated) <- "integrated"

# Save integrated seurat object
save(mixed.integrated, file="/Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza/integrate_mixed.RData")
# RunPCA and RunUMAP on integrated object-----------------------------------------------------------
mixed.integrated <- RunPCA(mixed.integrated, verbose = FALSE)
#Assess dimensionality by observing Elbow Plot and PCA Plot
ElbowPlot(mixed.integrated)
PCAPlot(mixed.integrated,group.by = "sample")
# Run UMAP, FindNeighbors, FindClusters on integrated object for 15 reductions
mixed.integrated <- RunUMAP(mixed.integrated, reduction = "pca", dims = 1:15)
mixed.integrated <- FindNeighbors(mixed.integrated, dims = 1:15)
mixed.integrated <- FindClusters(mixed.integrated,dims = 1:15,
                                resolution = c(0.2,0.3,0.6))
DimPlot(mixed.integrated, reduction = "umap") + ggtitle("Integrated Co-Cultured Datasets")
DimPlot(mixed.integrated, reduction = "umap", group.by = "CellLine.Treatment",
        cols = c('SNU719_Aza' = 'tomato2','SNU719_DMSO' = 'lightpink1','NCC24_Aza' = 'forestgreen', 'NCC24_DMSO' = 'seagreen3',
                 '1_Aza' = 'gray70', '1_DMSO' = 'gray70', '2_Aza' = 'gray70','2_DMSO' = 'gray70', '4_Aza' = 'gray70', '4_DMSO' = 'gray70', '5_Aza' = 'gray70', '5_DMSO' = 'gray70','6_Aza' = 'gray70', '6_DMSO' = 'gray70','7_Aza' = 'gray70', '7_DMSO' = 'gray70')) +
  ggtitle("Integrated Co-Cultured Datasets")
DimPlot(mixed.integrated, group.by = "CellLine", split.by = "Treatment", ncol = 2)
DimPlot(mixed.integrated, group.by = "CellLine.Treatment")
# Save integrated seurat object
save(mixed.integrated, file="/Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza/integrate_mixed_withUMAP_withIdents_23March.RData")

# Normalize data in RNA assay ----------------------------------------------------------
DefaultAssay(mixed.integrated) <- "RNA"
Idents(mixed.integrated)<-mixed.integrated$CellLine.Treatment
mixed.integrated <- NormalizeData(mixed.integrated, verbose = FALSE)
mixed.integrated <- FindVariableFeatures(mixed.integrated, selection.method = "vst", nfeatures = "3000")
top10 <- head(VariableFeatures(mixed.integrated), 10)
p3 <- VariableFeaturePlot(mixed.integrated)
LabelPoints(plot = p3, points = top10, repel = TRUE)
mixed.integrated <- ScaleData(mixed.integrated, features = VariableFeatures(mixed.integrated), vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))

### Set Idents -----------------------------------------------------------
# Name clusters according to Cell Line
table(Idents(mixed.integrated))
new.cluster.ids <- c("1", "2", "4", "5", "6", "7", "NCC24", "NCC24", "SNU719","SNU719")
names(new.cluster.ids) <- levels(mixed.integrated)
mixed.integrated <- RenameIdents(mixed.integrated, new.cluster.ids)
# Update CellLine
mixed.integrated$CellLine <- Idents(mixed.integrated)
# CellLine.Treatment
Idents(mixed.integrated)<- mixed.integrated$CellLine 
mixed.integrated$CellLine.Treatment<- paste(Idents(mixed.integrated), mixed.integrated$Treatment, sep = "_")
# CellLine.Culture
Idents(mixed.integrated)<-mixed.integrated$CellLine
mixed.integrated$CellLine.Culture <- paste(Idents(mixed.integrated), mixed.integrated$Culture, sep = "_")
# CellLine.Treatment.Culture
Idents(mixed.integrated)<-mixed.integrated$CellLine.Treatment
mixed.integrated$CellLine.Treatment.Culture <- paste(Idents(mixed.integrated), mixed.integrated$Culture, sep = "_")

DefaultAssay(mixed.integrated) <- "RNA"
Idents(mixed.integrated)<-mixed.integrated$CellLine.Treatment
# Save normalized, scaled object with corrected idents
save(mixed.integrated, file="/Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza/integrate_mixed_scaledRNAassay.RData")

# Average Expression-----------------------------------------------------------
Idents(mixed.integrated) <- mixed.integrated$CellLine
SNU719.cells <- subset(mixed.integrated, idents = c("SNU719"))
Idents(SNU719.cells) <- mixed.integrated$Treatment
avg.SNU719.cells <- log1p(AverageExpression(SNU719.cells, verbose = FALSE)$RNA)
p1 <- ggplot(avg.SNU719.cells, aes(`DMSO`, `Aza`)) + geom_point() + ggtitle("SNU719-DMSO vs SNU719-Aza in Co-Culture")
HoverLocator(plot = `p1`)
# genes.to.label = c("LGALS1", "IFI27", "MIF", "ADIRF", "HOXA10", "MAP1LC3A","ID4","TPD52L1")
# LabelPoints(p1, points = genes.to.label)

NCC24.cells <- subset(mixed.integrated, idents = c("NCC24"))
Idents(NCC24.cells) <- "Treatment"
avg.NCC24.cells <- log1p(AverageExpression(NCC24.cells, verbose = FALSE)$RNA)
p2 <- ggplot(avg.NCC24.cells, aes(`DMSO`, `Aza`)) + geom_point() + ggtitle("NCC24-DMSO vs NCC24-Aza in Co-Culture")
HoverLocator(plot = `p2`)
# genes.to.label = c("DKK1","PAGE4","AREG", "CAMK2N1","NFIA", "CCL20", "CXCL10", "TFF1", "CTNNB1","LCN1", "HIST1H4C","SDC4")
# LabelPoints(p2, points = genes.to.label)

#save objects of SNU and NCC 
save(SNU719.cells, file = "/Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza/co-SNU719-cells.Rdata")
save(NCC24.cells, file = "/Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza/co-NCC24-cells.Rdata")



# Identify conserved cell line markers  ----------------------------------------------------------
DefaultAssay(mixed.integrated) <- "RNA"
Idents(mixed.integrated) <- "CellLine"
# Identify genes that are conserved markers irrespective of treatment with Aza in the SNU719 and NCC24 clusters
SNU719.markers <- FindConservedMarkers(mixed.integrated, ident.1 = 'SNU719', grouping.var = "Treatment")
head(SNU719.markers)
NCC24.markers <- FindConservedMarkers(mixed.integrated, ident.1 = 'NCC24', grouping.var = "Treatment")
head(NCC24.markers)
# Explore the marker genes for the SNU719 and NCC24 cluster in FP
# Marker genes: High avglog_FC, pct1 >> pct2.
FeaturePlot(mixed.integrated, features = c("ALDH3A1", "PARK7", "LYZ", "S100P"))

Idents(mixed.integrated) <- "CellLine"
markers.to.plot <- c("ALDH3A1","AGR2", "PARK7", "S100P", "SPINK1","ODAM")
DotPlot(mixed.integrated, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "Treatment") + RotatedAxis()


# DE Analaysis: Comparing DMSO and Aza for SNU719 and NCC24 each ------------------------------------------------------------------
# Load annotations
annotations <- read.csv("/users/serena/scrnaseq/data/annotation.csv")
# Create new object to perform analysis on:
mixed.integrated.new <- mixed.integrated
Idents(mixed.integrated.new)<-mixed.integrated.new$CellLine.Treatment
DefaultAssay(mixed.integrated.new) <- "RNA"

### SNU719
# FindMarkers
mixed.SNU719.Treatment.markers  <- FindMarkers(mixed.integrated.new, ident.1 = "SNU719_Aza", ident.2 = "SNU719_DMSO", verbose = FALSE)
mixed.SNU719.Treatment.markers <- mixed.SNU719.Treatment.markers %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

# save list of de genes
write.csv(mixed.SNU719.Treatment.markers,file = "/Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza/All_DE_genes_SNU719_mixedDMSOvmixedAza.csv")
save(mixed.SNU719.Treatment.markers, file = "/Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza/All_DE_genes_SNU719_mixedDMSOvmixedAza.Rdata")

### NCC24
# FindMarkers
mixed.NCC24.Treatment.markers  <- FindMarkers(mixed.integrated.new, ident.1 = "NCC24_Aza", ident.2 = "NCC24_DMSO", verbose = FALSE)
mixed.NCC24.Treatment.markers <- mixed.NCC24.Treatment.markers %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

# save list of de genes
write.csv(mixed.NCC24.Treatment.markers,file = "/Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza/All_DE_genes_NCC24_mixedDMSOvmixedAza.csv")
save(mixed.NCC24.Treatment.markers, file = "/Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza/All_DE_genes_NCC24_mixedDMSOvmixedAza.Rdata")


# Visualizing DE Genes  -----------------------------------------------------------
Idents(mixed.integrated.new)<-mixed.integrated.new$CellLine.Treatment
DefaultAssay(mixed.integrated.new) <- "RNA"

### SNU719
# subset top 20 DE genes for each treatment condition
top20.mixed.SNU719.Treatment.markers  <- (top_n(mixed.SNU719.Treatment.markers, n = 20, wt = avg_logFC))$gene
top20.mixed.SNU719.Treatment.markers <- append(top20.mixed.SNU719.Treatment.markers, (top_n(mixed.SNU719.Treatment.markers, n = 20, wt = -avg_logFC))$gene)
# Create object of just SNU719_DMSO and SNU719_Aza cells
mixed.SNU719 <- subset(mixed.integrated.new, idents = c("SNU719_DMSO", "SNU719_Aza"))
# plot heatmap of top 20 for each treatment
DoHeatmap(mixed.SNU719, group.by = "CellLine.Treatment",features = top20.mixed.SNU719.Treatment.markers, 
          group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
  scale_fill_distiller(palette = "RdBu")

### NCC24
# subset top 20 DE genes for each treatment condition
top20.mixed.NCC24.Treatment.markers  <- (top_n(mixed.NCC24.Treatment.markers, n = 20, wt = avg_logFC))
top20.mixed.NCC24.Treatment.markers <- append(top20.mixed.NCC24.Treatment.markers, (top_n(mixed.NCC24.Treatment.markers, n = 20, wt = -avg_logFC))$gene)
# Create object of just NCC24_DMSO and NCC24_Aza cells
mixed.NCC24 <- subset(mixed.integrated.new, idents = c("NCC24_DMSO", "NCC24_Aza"))
# plot heatmap of top 20 for each treatment
DoHeatmap(mixed.NCC24, group.by = "CellLine.Treatment",features = top20.mixed.NCC24.Treatment.markers, 
          group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu")


# Feature Plot of selected genes across treatment condition
# modify genes under 'features'
# moddify object to be plotted as desiredzw
FeaturePlot(mixed.SNU719, features = c("IFI27","IFI6"," LGALS1","BMP4"), 
            split.by = "Treatment", max.cutoff = 3, cols = c("grey", "red"))

# Violin plot of selected genes across treatment condition
# modify genes under 'features'
plots <- VlnPlot(mixed.SNU719, 
                 features = c("IFI27","IFI6","LGALS1", "RPL17"), 
                 split.by = "Treatment", 
                 group.by = 'CellLine',
                 pt.size = 0, 
                 combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x +
                  theme(axis.title = element_text(size=8),
                        axis.text.x = element_text(size=8,angle = 0, hjust = 0.5),
                        axis.text.y = element_text(size=8),
                        plot.title=element_text(size=8,face="bold"),
                        legend.position = 'top'))
CombinePlots(plots = plots, nrow = 2)


# # Old Code -----------------------------------------------------------
# 
# ### SNU719 cells: DMSO v Aza
# # Switch current ident to the column holding Treatment Label
# Idents(SNU719.cells)<-SNU719.cells$Treatment
# # Set assay to either RNA or integrated
# DefaultAssay(SNU719.cells) <- "integrated"
# 
# # Find Marker Genes between SNU719-DMSO and SNU719-Aza
# co.SNU719.Treatment.markers <- FindMarkers(SNU719.cells, ident.1= 'Aza', ident.2 = "DMSO")
# # co.SNU719.Treatment.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# 
#  # Combine markers with gene descriptions
# co.SNU719.Treatment.markers_ann <- co.SNU719.Treatment.markers %>%
#   rownames_to_column(var="gene") %>%
#   left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
# 
# # Save before continuing any further
# write.csv(co.SNU719.Treatment.markers_ann,"Markers for co SNU719 DMSO vs co SNU719 Aza.csv")
# save(co.SNU719.Treatment.markers, file = "/Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza/Markers-coSNU719DMSO-v-coSNU719Aza.Rdata")
# save(co.SNU719.Treatment.markers_ann, file = "/Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza/Markers-coSNU719DMSO-v-coSNU719Aza_ann.Rdata")
# 
# ### NCC24 cells: DMSO v Aza
# # Switch current ident to the column holding Treatment Label
# Idents(NCC24.cells)<-NCC24.cells$Treatment
# 
# # Find Marker Genes between NCC24-DMSO and NCC24-Aza
# co.NCC24.Treatment.markers <- FindMarkers(NCC24.cells, ident.1= 'Aza', ident.2 = "DMSO")
# # co.NCC24.Treatment.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# 
# # Combine markers with gene descriptions
# co.NCC24.Treatment.markers_ann <- co.NCC24.Treatment.markers %>%
#   rownames_to_column(var="gene") %>%
#   left_join(y = unique(annotations[, c("gene_name", "description")]),
#             by = c("gene" = "gene_name"))
# 
# # Save before continuing any further
# write.csv(co.NCC24.Treatment.markers_ann,"Markers for co NCC24 DMSO vs co NCC24 Aza.csv")
# save(co.NCC24.Treatment.markers, file = "/Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza/Markers-coNCC24DMSO-v-coNCC24Aza.Rdata")
# save(co.NCC24.Treatment.markers_ann, file = "/Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza/Markers-coNCC24DMSO-v-coNCC24Aza_ann.Rdata")
# De:Avg of Each Cluster
# # replace spaces with underscores '_' so ggplot2 doesn't fail
# orig.levels <- levels(mixed.SNU719)
# Idents(mixed.SNU719) <- gsub(pattern = " ", replacement = "_", x = Idents(mixed.SNU719))
# orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
# levels(mixed.SNU719) <- orig.levels
# cluster.averages <- AverageExpression(mixed.SNU719, return.seurat = TRUE)
# DoHeatmap(cluster.averages, features = genes.SNU719.Aza, group.bar = TRUE, draw.lines = FALSE, size = 4, hjust = 1, angle = 0) +
#   scale_fill_distiller(palette = "RdBu")
#Avg of Each Cluster
# replace spaces with underscores '_' so ggplot2 doesn't fail
# orig.levels <- levels(mixed.NCC24)
# Idents(mixed.NCC24) <- gsub(pattern = " ", replacement = "_", x = Idents(mixed.NCC24))
# orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
# levels(mixed.NCC24) <- orig.levels
# cluster.averages <- AverageExpression(mixed.NCC24, return.seurat = TRUE)
# DoHeatmap(cluster.averages, features = genes.NCC24.Aza, group.bar = TRUE, draw.lines = FALSE, size = 4, hjust = 1, angle = 0) +
#   scale_fill_distiller(palette = "RdBu")

