#!/usr/bin/env Rscript

# README ------------------------------------------------------------------

# Script version: 0
# Question: How do culture conditions affect the expression of DMSO cells?

# Set-up packages, working directory, datasets, set idents -----------------------------------------------------------
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
setwd("/Users/serena/scrnaseq/results/analysis/pureDMSOvcoDMSO")

## Load the datasets: 2 SNU719-DMSO, 1 NCC24-DMSO,1 mixed_DMSO_labeled
load("/Users/serena/scrnaseq/results/pure_indiv/april2017/SNU719-DMSO-april2017/v1/SCT-SNU719-DMSO-april2017-v1.RData")
load("/Users/serena/scrnaseq/results/pure_indiv/july2017/SNU719-DMSO-july2017/v1/SCT-SNU719-DMSO-july2017-v1.RData")
load("/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/mixed_DMSO_labeled_v0.RData")
load("/Users/serena/scrnaseq/results/pure_indiv/sept2017/NCC24-DMSO-sept2017/v1/SCT-NCC24-DMSO-sept2017-v1.RData")

### Set Idents for SNU719.DMSO.april2017.seuratobject
Idents(SNU719.DMSO.april2017.seuratobject)<-SNU719.DMSO.april2017.seuratobject$CellLine
# CellLine.Treatment
SNU719.DMSO.april2017.seuratobject$CellLine.Treatment<- paste(Idents(SNU719.DMSO.april2017.seuratobject), SNU719.DMSO.april2017.seuratobject$Treatment, sep = "_")
# CellLine.Culture, CellLine.Treatment.Culture
SNU719.DMSO.april2017.seuratobject$CellLine.Culture <- "SNU719_Mono"
SNU719.DMSO.april2017.seuratobject$CellLine.Treatment.Culture <- "SNU719_DMSO_Mono"
# Set Idents
Idents(SNU719.DMSO.april2017.seuratobject) <-SNU719.DMSO.april2017.seuratobject$CellLine.Culture
# save
save(SNU719.DMSO.april2017.seuratobject, file = "/Users/serena/scrnaseq/results/analysis/pureDMSOvcoDMSO/SCT_SNU719_DMSO_april2017_correct_Idents.RData")

### Set Idents for SNU719.DMSO.july2017.seuratobject
Idents(SNU719.DMSO.july2017.seuratobject)<-SNU719.DMSO.july2017.seuratobject$CellLine
# Add in Idents of CellLine.Treatment
SNU719.DMSO.july2017.seuratobject$CellLine.Treatment<- paste(Idents(SNU719.DMSO.july2017.seuratobject), SNU719.DMSO.july2017.seuratobject$Treatment, sep = "_")
# Add in Idents of CellLine.Culture
SNU719.DMSO.july2017.seuratobject$CellLine.Culture <- "SNU719_Mono"
SNU719.DMSO.july2017.seuratobject$CellLine.Treatment.Culture <- "SNU719_DMSO_Mono"
# Set Idents
Idents(SNU719.DMSO.july2017.seuratobject) <- SNU719.DMSO.july2017.seuratobject$CellLine.Culture
# save
save(SNU719.DMSO.july2017.seuratobject, file = "/Users/serena/scrnaseq/results/analysis/pureDMSOvcoDMSO/SCT_SNU719_DMSO_july2017_correct_Idents.RData")

### Set Idents for NCC24.DMSO.sept2017.seuratobject
Idents(NCC24.DMSO.sept2017.seuratobject)<-NCC24.DMSO.sept2017.seuratobject$CellLine
# Add in Idents of CellLine.Treatment
NCC24.DMSO.sept2017.seuratobject$CellLine.Treatment<- paste(Idents(NCC24.DMSO.sept2017.seuratobject), NCC24.DMSO.sept2017.seuratobject$Treatment, sep = "_")
# Add in Idents of CellLine.Culture
NCC24.DMSO.sept2017.seuratobject$CellLine.Culture <- "NCC24_Mono"
NCC24.DMSO.sept2017.seuratobject$CellLine.Treatment.Culture <- "NCC24_DMSO_Mono"
#Set Idents
Idents(NCC24.DMSO.sept2017.seuratobject) <- NCC24.DMSO.sept2017.seuratobject$CellLine.Culture
# save
save(NCC24.DMSO.sept2017.seuratobject, file = "/Users/serena/scrnaseq/results/analysis/pureDMSOvcoDMSO/SCT_NCC24_DMSO_sept2017_correct_Idents.RData")

### Set Idents for mixed_DMSO_labeled
# Check names of clusters
table(Idents(mixed_DMSO_labeled))
# Name clusters according to Cell Line
new.cluster.ids <- c("SNU719", "1", "2", "NCC24", "4", "5", "6", "7")
names(new.cluster.ids) <- levels(mixed_DMSO_labeled)
mixed_DMSO_labeled <- RenameIdents(mixed_DMSO_labeled, new.cluster.ids)
# Change CellLine metadata in idents
mixed_DMSO_labeled$CellLine <- Idents(mixed_DMSO_labeled)
# Add in Idents of CellLine.Treatment
mixed_DMSO_labeled$CellLine.Treatment<- paste(Idents(mixed_DMSO_labeled), mixed_DMSO_labeled$Treatment, sep = "_")
# Add in Idents of CellLine.Culture
mixed_DMSO_labeled$CellLine.Culture<- paste(Idents(mixed_DMSO_labeled), mixed_DMSO_labeled$Culture, sep = "_")
Idents(mixed_DMSO_labeled) <- mixed_DMSO_labeled$CellLine.Treatment
mixed_DMSO_labeled$CellLine.Treatment.Culture <- paste(Idents(mixed_DMSO_labeled), mixed_DMSO_labeled$Culture, sep = "_")
# Save correct Idents mixed_DMSO_labeled
save(mixed_DMSO_labeled, file = "/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/mixed_DMSO_labeled_correct_Idents.RData")
save(mixed_DMSO_labeled, file = "/Users/serena/scrnaseq/results/analysis/pureDMSOvcoDMSO/mixed_DMSO_labeled_correct_Idents.RData")

# Integrate datasets: SNU719 and mixed_DMSO--------------------------------------
mixed.SNU.DMSO.list <- list(SNU719.DMSO.april2017.seuratobject,SNU719.DMSO.july2017.seuratobject, mixed_DMSO_labeled)
mixed.SNU.DMSO.features <- SelectIntegrationFeatures(object.list = mixed.SNU.DMSO.list, nfeatures = 3000)
mixed.SNU.DMSO.list <- PrepSCTIntegration(object.list = mixed.SNU.DMSO.list, anchor.features = mixed.SNU.DMSO.features)
mixed.SNU.DMSO.anchors <- FindIntegrationAnchors(object.list = mixed.SNU.DMSO.list, normalization.method = "SCT", anchor.features = mixed.SNU.DMSO.features)
mixed.SNU.DMSO.integrated <- IntegrateData(anchorset = mixed.SNU.DMSO.anchors, normalization.method = "SCT")
DefaultAssay(mixed.SNU.DMSO.integrated) <- "integrated"

# save 
save(mixed.SNU.DMSO.integrated, file = "/Users/serena/scrnaseq/results/analysis/pureDMSOvcoDMSO/mixed_SNU_DMSO_integrated.Rdata")

# RunPCA Elbow Plot UMAP
mixed.SNU.DMSO.integrated <- RunPCA(mixed.SNU.DMSO.integrated, verbose = FALSE)
ElbowPlot(mixed.SNU.DMSO.integrated)
PCAPlot(mixed.SNU.DMSO.integrated,group.by = "sample")
mixed.SNU.DMSO.integrated <- RunUMAP(mixed.SNU.DMSO.integrated, reduction = "pca", dims = 1:15)

DimPlot(mixed.SNU.DMSO.integrated, reduction = "umap") + ggtitle("Integrated Mono and Co SNU719 DMSO Datasets")

DimPlot(mixed.SNU.DMSO.integrated, group.by = "CellLine.Culture",
        cols = c('SNU719_Co' = 'tomato2', 'SNU719_Mono' = 'lightpink1',
                 'NCC24_Co' = "gray70",'1_Co' = "gray70",'2_Co' = "gray70",'4_Co' = "gray70",'5_Co' = "gray70",'6_Co' = "gray70",'7_Co' = "gray70")) +
  ggtitle("Integrated Mono and Co SNU719 DMSO Datasets")

DimPlot(mixed.SNU.DMSO.integrated, group.by = "CellLine",
        cols = c('SNU719' = 'tomato2',
                 'NCC24' = "forestgreen",'1' = "gray70",'2' = "gray70",'4' = "gray70",'5' = "gray70",'6' = "gray70",'7' = "gray70")) +
  ggtitle("Integrated Mono and Co SNU719 DMSO Datasets")

Idents(mixed.SNU.DMSO.integrated) <- "CellLine"
DimPlot(mixed.SNU.DMSO.integrated, split.by = "Culture")
        

# Save integrated seurat object
save(mixed.SNU.DMSO.integrated, file="/Users/serena/scrnaseq/results/analysis/pureDMSOvcoDMSO/integrate_mixed_SNU_DMSO_withUMAP.RData")

# Set-up SNU object for DE: Normalize, Set Idents, ScaleData  ----------------------------------------------------------
# Set assay to RNA
DefaultAssay(mixed.SNU.DMSO.integrated) <- "RNA"
# Normalize RNA data for visualization purposes
mixed.SNU.DMSO.integrated <- NormalizeData(mixed.SNU.DMSO.integrated, verbose = FALSE)
# Scale Data, and save the object again
# FVF to plot 
mixed.SNU.DMSO.integrated <- FindVariableFeatures(mixed.SNU.DMSO.integrated, selection.method = "vst", nfeatures = "3000")
top10 <- head(VariableFeatures(mixed.SNU.DMSO.integrated), 10)
p3 <- VariableFeaturePlot(mixed.SNU.DMSO.integrated)
LabelPoints(plot = p3, points = top10, repel = TRUE)
mixed.SNU.DMSO.integrated <- ScaleData(mixed.SNU.DMSO.integrated, 
                                       features = VariableFeatures(mixed.SNU.DMSO.integrated), 
                                       vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
# Save normalized, scaled object
save(mixed.SNU.DMSO.integrated, file="/Users/serena/scrnaseq/results/analysis/pureDMSOvcoDMSO/integrate_mixed_SNU_DMSO_withUMAP_and_ScaleData.RData")

# FCM: SNU Across Culture--------------------------------------
DefaultAssay(SNU719.Culture.Cells) <- "RNA"
Idents(SNU719.Culture.Cells) <- "CellLine"
# #Find conserved genes for SNU719-DMSO across culture (Co v Mono)
SNU719.markers <- FindConservedMarkers(SNU719.Culture.Cells, ident.1 = 'SNU719', grouping.var = "Culture", verbose = FALSE)
head(SNU719.markers)

# Explore the marker genes for  SNU719 in FP
# Marker genes: High avglog_FC, pct1 >> pct2.
FeaturePlot(mixed.SNU.DMSO.integrated, features = c("SPINK4"))

Idents(mixed.SNU.DMSO.integrated) <- "CellLine"
markers.to.plot <- c("HIST1H4C","CCK","MT2A","LYZ")
DotPlot(mixed.SNU.DMSO.integrated, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "Culture") + RotatedAxis()


# Average Expression of SNU719-----------------------------------------------------------
Idents(mixed.SNU.DMSO.integrated) <- mixed.SNU.DMSO.integrated$CellLine
SNU719.cells <- subset(mixed.SNU.DMSO.integrated, idents = c("SNU719"))
Idents(SNU719.cells) <-mixed.SNU.DMSO.integrated$Culture
avg.SNU719.cells <- log1p(AverageExpression(SNU719.cells, verbose = FALSE)$RNA)
p1 <- ggplot(avg.SNU719.cells, aes(`Mono`, `Co`)) + geom_point() + ggtitle("DMSO-treated SNU719 in Mono vs DMSO-treated SNU719 in Co-Culture")
HoverLocator(plot = `p1`)
genes.to.label = c("IFI6", "IFI27", "STAT1", "NOP53", "SLC25A6", "PRDX11","ACB","TPD52L1")
LabelPoints(p1, points = genes.to.label)
#DE: SNU--------------------------------------
# Load annotations
annotations <- read.csv("/users/serena/scrnaseq/data/annotation.csv")
# Create new object to perform analysis on:
mixed.SNU.DMSO.integrated.new <- mixed.SNU.DMSO.integrated
Idents(mixed.SNU.DMSO.integrated.new)<-mixed.SNU.DMSO.integrated.new$CellLine.Culture
DefaultAssay(mixed.SNU.DMSO.integrated.new) <- "RNA"

# FindMarkers
DMSO.SNU719.Culture.response <- FindMarkers(mixed.SNU.DMSO.integrated.new, ident.1 = "SNU719_Co", ident.2 = "SNU719_Mono", verbose = TRUE, min.diff.pct = 0.3)
DMSO.SNU719.Culture.response <- DMSO.SNU719.Culture.response %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

# save list of de genes
write.csv(DMSO.SNU719.Culture.response,file = "/Users/serena/scrnaseq/results/analysis/pureDMSOvcoDMSO/All_DE_genes_SNU719_pureDMSOvcoDMSO.csv")
save(DMSO.SNU719.Culture.response, file = "/Users/serena/scrnaseq/results/analysis/pureDMSOvcoDMSO/All_DE_genes_SNU719_pureDMSOvcoDMSO.Rdata")

# Visualizing DE Genes  -----------------------------------------------------------
# List top 20 up reg and top 20 down reg
top20.DMSO.SNU719.Culture.markers <- (top_n(DMSO.SNU719.Culture.response, n = 20, wt = avg_logFC))$gene
top20.DMSO.SNU719.Culture.markers <- append(top20.DMSO.SNU719.Culture.markers, (top_n(DMSO.SNU719.Culture.response, n = 20, wt = -avg_logFC))$gene)

# Create object of just SNU719_Mono and SNU719_Culture cells
DMSO.SNU719.Cells <- subset(mixed.SNU.DMSO.integrated.new, idents = c("SNU719_Mono", "SNU719_Co"))

# Plot Heatmap of these genes
DefaultAssay(DMSO.SNU719.Cells) <- "RNA"
DoHeatmap(DMSO.SNU719.Cells, group.by = "CellLine.Culture",features = top20.DMSO.SNU719.Culture.markers, 
          group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
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
                  theme(axis.title = element_text(size=8),
                        axis.text.x = element_text(size=8,angle = 0, hjust = 0.5),
                        axis.text.y = element_text(size=8),
                        plot.title=element_text(size=8,face="bold"),
                        legend.position = 'top'))
CombinePlots(plots = plots, nrow = 2)
# OLD CODE  ----------------------------------------------------------
# Integrate datasets: NCC24 and mixed_DMSO--------------------------------------
Idents(NCC24.DMSO.sept2017.seuratobject) <- "CellLine.Culture"
Idents(mixed_DMSO_labeled) <- "CellLine.Culture"
mixed.NCC.DMSO.list <- list(NCC24.DMSO.sept2017.seuratobject, mixed_DMSO_labeled)
mixed.NCC.DMSO.features <- SelectIntegrationFeatures(object.list = mixed.NCC.DMSO.list, nfeatures = 3000)
mixed.NCC.DMSO.list <- PrepSCTIntegration(object.list = mixed.NCC.DMSO.list, anchor.features = mixed.NCC.DMSO.features)
mixed.NCC.DMSO.anchors <- FindIntegrationAnchors(object.list = mixed.NCC.DMSO.list, normalization.method = "SCT", anchor.features = mixed.NCC.DMSO.features)
mixed.NCC.DMSO.integrated <- IntegrateData(anchorset = mixed.NCC.DMSO.anchors, normalization.method = "SCT")

# Set assay
DefaultAssay(mixed.NCC.DMSO.integrated) <- "integrated"

# save in case 4.07 pm
save(mixed.NCC.DMSO.integrated, file = "mixed_NCC_DMSO_integrated.Rdata")

# RunPCA Elbow Plot UMAP
mixed.NCC.DMSO.integrated <- RunPCA(mixed.NCC.DMSO.integrated, verbose = FALSE)
png("Elbow Plot Integrating mixedDMSO and pureNCC-DMSO datasets.png")
ElbowPlot(mixed.NCC.DMSO.integrated)
dev.off()
png("PCA Plot Integrating mixedDMSO and pureNCC-DMSO datasets.png")
PCAPlot(mixed.NCC.DMSO.integrated,group.by = "sample")
dev.off()
mixed.NCC.DMSO.integrated <- RunUMAP(mixed.NCC.DMSO.integrated, reduction = "pca", dims = 1:15)

table(Idents(mixed.NCC.DMSO.integrated))

DimPlot(mixed.NCC.DMSO.integrated, reduction = "umap") + ggtitle("Integrated Mono and Co NCC24 DMSO Datasets")
DimPlot(mixed.NCC.DMSO.integrated, group.by = "CellLine.Culture",
        cols = c('NCC24_Co' = 'forestgreen', 'NCC24_Mono' = 'seagreen3',
                 'SNU719_Co' = "tomato2",'1_Co' = "gray70",'2_Co' = "gray70",'4_Co' = "gray70",'5_Co' = "gray70",'6_Co' = "gray70",'7_Co' = "gray70"))

DimPlot(mixed.NCC.DMSO.integrated, group.by = "CellLine.Culture",
        cols = c('NCC24_Co' = 'forestgreen', 'NCC24_Mono' = 'seagreen3',
                 'SNU719_Co' = "gray70",'1_Co' = "gray70",'2_Co' = "gray70",'4_Co' = "gray70",'5_Co' = "gray70",'6_Co' = "gray70",'7_Co' = "gray70"))


DimPlot(mixed.NCC.DMSO.integrated, group.by = "CellLine",
        cols = c('SNU719' = 'tomato2',
                 'NCC24' = "forestgreen",'1' = "gray70",'2' = "gray70",'4' = "gray70",'5' = "gray70",'6' = "gray70",'7' = "gray70")) +
  ggtitle("Integrated Mono and Co SNU719 DMSO Datasets")

# Save integrated seurat object
save(mixed.NCC.DMSO.integrated, file="/Users/serena/scrnaseq/results/analysis/pureDMSOvcoDMSO/integrate_mixed_NCC_DMSO_withUMAP.RData")


# Set-up NCC object for DE: Normalize, Set Idents, ScaleData  ----------------------------------------------------------
# Set assay to RNA
DefaultAssay(mixed.NCC.DMSO.integrated) <- "RNA"
# Normalize RNA data for visualization purposes
mixed.NCC.DMSO.integrated <- NormalizeData(mixed.NCC.DMSO.integrated, verbose = FALSE)
# Scale Data, and save the object again
# FVF to plot 
mixed.NCC.DMSO.integrated <- FindVariableFeatures(mixed.NCC.DMSO.integrated, selection.method = "vst", nfeatures = "3000")
top10 <- head(VariableFeatures(mixed.NCC.DMSO.integrated), 10)
p3 <- VariableFeaturePlot(mixed.NCC.DMSO.integrated)
LabelPoints(plot = p3, points = top10, repel = TRUE)

# ScaleData for Heatmap
mixed.NCC.DMSO.integrated <- ScaleData(mixed.NCC.DMSO.integrated, vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))

# Save normalized, scaled object
save(mixed.NCC.DMSO.integrated, file="/Users/serena/scrnaseq/results/analysis/pureDMSOvcoDMSO/integrate_mixed_NCC_DMSO_withUMAP_and_ScaleData.RData")

# FCM: NCC Across Culture--------------------------------------
DefaultAssay(mixed.NCC.DMSO.integrated) <- "RNA"
Idents(mixed.NCC.DMSO.integrated) <- "CellLine.Culture"
# #Find conserved genes for SNU719-DMSO across culture (Co v Mono)
NCC24.markers <- FindConservedMarkers(mixed.NCC.DMSO.integrated, ident.1 = 'NCC24_Co', grouping.var = "Treatment")
head(NCC24.markers)

# Explore the marker genes for  SNU719 in FP
# Marker genes: High avglog_FC, pct1 >> pct2.
FeaturePlot(mixed.NCC.DMSO.integrated, features = c())

Idents(mixed.NCC.DMSO.integrated) <- "CellLine"
markers.to.plot <- c()
DotPlot(mixed.NCC.DMSO.integrated, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "Culture") + RotatedAxis()
# DE: NCC--------------------------------------
# Load annotations
annotations <- read.csv("/users/serena/scrnaseq/data/annotation.csv")
# Create new object to perform analysis on:
mixed.NCC.DMSO.integrated.new <- mixed.NCC.DMSO.integrated
# Set Idents - What idents to set to?
# I am testing NCC24_Co vs NCC24_Mono
Idents(mixed.NCC.DMSO.integrated.new)<-mixed.NCC.DMSO.integrated.new$CellLine.Culture
# Set RNA assay
DefaultAssay(mixed.NCC.DMSO.integrated.new) <- "RNA"

# FindMarkers
NCC24.Culture.response <- FindMarkers(mixed.NCC.DMSO.integrated.new, ident.1 = "NCC24_Co", ident.2 = "NCC24_Mono", verbose = TRUE)
NCC24.Culture.response <- NCC24.Culture.response %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

# List top 20 up reg and top 20 down reg
genes.NCC24.Culture <- (top_n(NCC24.Culture.response, n = 20, wt = avg_logFC))$gene
genes.NCC24.Culture <- append(genes.NCC24.Culture, (top_n(NCC24.Culture.response, n = 20, wt = -avg_logFC))$gene)

# Create object of just NCC24_DMSO and NCC24_Culture cells
NCC24.Culture.Cells <- subset(mixed.NCC.DMSO.integrated.new, idents = c("NCC24_Mono", "NCC24_Co"))

# Plot Heatmap of these genes
DefaultAssay(NCC24.Culture.Cells) <- "RNA"
DoHeatmap(NCC24.Culture.Cells, group.by = "CellLine.Culture",features = genes.NCC24.Culture, group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu")
# Perform FCM v2: Integrated mixed, NCC, SNU datasets--------------------------------------
Idents(NCC24.DMSO.sept2017.seuratobject) <- "CellLine.Culture"
Idents(mixed_DMSO_labeled) <- "CellLine.Culture"
Idents(SNU719.DMSO.april2017.seuratobject) <- "CellLine.Culture"
Idents(SNU719.DMSO.july2017.seuratobject) <- "CellLine.Culture"
mixed.DMSO.list <- list(SNU719.DMSO.april2017.seuratobject,SNU719.DMSO.july2017.seuratobject, NCC24.DMSO.sept2017.seuratobject, mixed_DMSO_labeled)
mixed.DMSO.features <- SelectIntegrationFeatures(object.list = mixed.DMSO.list, nfeatures = 3000)
mixed.DMSO.list <- PrepSCTIntegration(object.list = mixed.DMSO.list, anchor.features = mixed.DMSO.features)
mixed.DMSO.anchors <- FindIntegrationAnchors(object.list = mixed.DMSO.list, normalization.method = "SCT", anchor.features = mixed.DMSO.features)
mixed.DMSO.integrated <- IntegrateData(anchorset = mixed.DMSO.anchors, normalization.method = "SCT")

DefaultAssay(mixed.DMSO.integrated) <- "RNA"
mixed.DMSO.integrated <- NormalizeData(mixed.DMSO.integrated, verbose = FALSE)
mixed.DMSO.integrated <- FindVariableFeatures(mixed.DMSO.integrated, nfeatures = 3000)
mixed.DMSO.integrated <- ScaleData(mixed.DMSO.integrated, features = VariableFeatures(mixed.DMSO.integrated), vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))

# Save normalized, scaled object
save(mixed.DMSO.integrated, file="/Users/serena/scrnaseq/results/analysis/pureDMSOvcoDMSO/integrate_mixed_SNU_NCC_DMSO_withUMAP_and_ScaleData.RData")

mixed.DMSO.integrated.new <- mixed.DMSO.integrated

DefaultAssay(mixed.DMSO.integrated.new) <- "RNA"
Idents(mixed.DMSO.integrated.new) <- "CellLine"

# Create object of just NCC24 and SNU719 cells
SNU719.NCC24.Culture.Cells <- subset(mixed.DMSO.integrated.new, idents = c("NCC24", "SNU719"))

# #Find conserved genes for SNU719-DMSO across culture (Co v Mono)
SNU719.conserved.markers <- FindConservedMarkers(mixed.DMSO.integrated.new, ident.1 = 'SNU719', grouping.var = "Culture", verbose = FALSE)
head(SNU719.conserved.markers)

NCC24.conserved.markers <- FindConservedMarkers(SNU719.NCC24.Culture.Cells, ident.1 = 'NCC24', grouping.var = "Culture", verbose = FALSE)
head(NCC24.conserved.markers)

# Explore conserved markers for SNU719 and NCC24
# Marker genes: High avglog_FC, pct1 >> pct2.
# Plot in all integrated
DefaultAssay(mixed.DMSO.integrated.new) <- "integrated"
mixed.DMSO.integrated.new <- RunPCA(mixed.DMSO.integrated.new, verbose = FALSE)
ElbowPlot(mixed.DMSO.integrated.new)
mixed.DMSO.integrated.new <- RunUMAP(mixed.DMSO.integrated.new, reduction = "pca", dims = 1:15)
DimPlot(mixed.DMSO.integrated.new, reduction = "umap")
FeaturePlot(mixed.DMSO.integrated.new, features = c("PHGR1", "ALDH3A1", "S100A14"))

Idents(mixed.DMSO.integrated.new) <- "CellLine"
markers.to.plot <- c("PHGR1", "ALDH3A1", "S100A14", "S100P", "TACC1", "ODAM")
DotPlot(mixed.DMSO.integrated.new, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "Culture")
DotPlot(mixed.DMSO.integrated.new, features = rev(markers.to.plot))

# Plot in subset of SNU719 and NCC24 cells only
DefaultAssay(mixed.DMSO.integrated.new) <- "RNA"
mixed.DMSO.integrated.new <- FindVariableFeatures(mixed.DMSO.integrated.new, selection.method = "vst", nfeatures = 3000)
mixed.DMSO.integrated.new <- ScaleData(mixed.DMSO.integrated.new)
mixed.DMSO.integrated.new <- RunPCA(mixed.DMSO.integrated.new, features = VariableFeatures(object = mixed.DMSO.integrated.new),  verbose = FALSE)
ElbowPlot(mixed.DMSO.integrated.new)
mixed.DMSO.integrated.new <- RunUMAP(mixed.DMSO.integrated.new, reduction = "pca", dims = 1:15)
FeaturePlot(SNU719.NCC24.Culture.Cells, features = c("CST1", "SPINK1", "LCN2"))

DefaultAssay(mixed.DMSO.integrated.new) <- "integrated"
mixed.DMSO.integrated.new <- RunPCA(mixed.DMSO.integrated.new,  verbose = FALSE)
ElbowPlot(mixed.DMSO.integrated.new)
mixed.DMSO.integrated.new <- RunUMAP(mixed.DMSO.integrated.new, reduction = "pca", dims = 1:15)

DimPlot(mixed.DMSO.integrated.new, group.by = "CellLine.Culture",
        cols = c('NCC24_Co' = 'forestgreen', 'NCC24_Mono' = 'seagreen3',
                 'SNU719_Co' = "tomato2",'SNU719_Mono'= 'lightpink1', '1_Co' = "gray70",'2_Co' = "gray70",'4_Co' = "gray70",'5_Co' = "gray70",'6_Co' = "gray70",'7_Co' = "gray70"))


# Explore conserved markers in mixed.SNU
Idents(mixed.SNU.DMSO.integrated) <- "CellLine"
markers.to.plot <- c("HIST1H4C","CCK","MT2A","LYZ")
DotPlot(mixed.SNU.DMSO.integrated, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "Culture") + RotatedAxis()

# Explore conserved markers in mixed.NCC
Idents(mixed.NCC.DMSO.integrated) <- "CellLine"
markers.to.plot <- c("HIST1H4C","CCK","MT2A","LYZ")
DotPlot(mixed.NCC.DMSO.integrated, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "Culture") + RotatedAxis()

# Explore UMAp
