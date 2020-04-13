#!/usr/bin/env Rscript

# README ------------------------------------------------------------------

# Script version: 0
# Question: How do culture conditions affect the expression of Aza cells?
# Only done for SNU719
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
setwd("/Users/serena/scrnaseq/results/analysis/pureAzavcoAza")

## Load the datasets: 2 SNU719-DMSO, 1 NCC24-DMSO,1 mixed_DMSO_labeled
load("/Users/serena/scrnaseq/results/pure_indiv/april2017/SNU719-Aza-april2017/v1/SCT-SNU719-Aza-april2017-v1.RData")
load("/Users/serena/scrnaseq/results/pure_indiv/sept2017/SNU719-Aza-sept2017/v1/SCT-SNU719-Aza-sept2017-v1.RData")
load("/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/mixed_Aza_labeled_v0.RData")

### Set Idents for SNU719.Aza.april2017.seuratobject
Idents(SNU719.Aza.april2017.seuratobject)<-SNU719.Aza.april2017.seuratobject$CellLine
# CellLine.Treatment
SNU719.Aza.april2017.seuratobject$CellLine.Treatment<- paste(Idents(SNU719.Aza.april2017.seuratobject), SNU719.Aza.april2017.seuratobject$Treatment, sep = "_")
# CellLine.Culture, CellLine.Treatment.Culture
SNU719.Aza.april2017.seuratobject$CellLine.Culture <- "SNU719_Mono"
SNU719.Aza.april2017.seuratobject$CellLine.Treatment.Culture <- "SNU719_Aza_Mono"
# Set Idents
Idents(SNU719.Aza.april2017.seuratobject) <- SNU719.Aza.april2017.seuratobject$CellLine.Culture
# save
save(SNU719.Aza.april2017.seuratobject, file = "/Users/serena/scrnaseq/results/analysis/pureAzavcoAza/SCT_SNU719_Aza_april2017_correct_Idents.RData")

### Set Idents for SNU719.Aza.sept2017.seuratobject
Idents(SNU719.Aza.sept2017.seuratobject)<-SNU719.Aza.sept2017.seuratobject$CellLine
# CellLine.Treatment
SNU719.Aza.sept2017.seuratobject$CellLine.Treatment<- paste(Idents(SNU719.Aza.sept2017.seuratobject), SNU719.Aza.sept2017.seuratobject$Treatment, sep = "_")
# CellLine.Culture, CellLine.Treatment.Culture
SNU719.Aza.sept2017.seuratobject$CellLine.Culture <- "SNU719_Mono"
SNU719.Aza.sept2017.seuratobject$CellLine.Treatment.Culture <- "SNU719_Aza_Mono"
# Set Idents
Idents(SNU719.Aza.sept2017.seuratobject) <- SNU719.Aza.sept2017.seuratobject$CellLine.Culture
# save
save(SNU719.Aza.sept2017.seuratobject, file = "/Users/serena/scrnaseq/results/analysis/pureAzavcoAza/SCT_SNU719_Aza_sept2017_correct_Idents.RData")

### Set Idents for mixed_Aza_labeled
# Check names of clusters
table(Idents(mixed_Aza_labeled))
# Name clusters according to Cell Line
new.cluster.ids <- c("SNU719", "1", "2", "NCC24", "4", "5", "6","7")
names(new.cluster.ids) <- levels(mixed_Aza_labeled)
mixed_Aza_labeled <- RenameIdents(mixed_Aza_labeled, new.cluster.ids)
# Change CellLine metadata in idents
mixed_Aza_labeled$CellLine <- Idents(mixed_Aza_labeled)
# Add in Idents of CellLine.Treatment
mixed_Aza_labeled$CellLine.Treatment<- paste(Idents(mixed_Aza_labeled), mixed_Aza_labeled$Treatment, sep = "_")
# Add in Idents of CellLine.Culture
mixed_Aza_labeled$CellLine.Culture<- paste(Idents(mixed_Aza_labeled), mixed_Aza_labeled$Culture, sep = "_")
Idents(mixed_Aza_labeled) <- mixed_Aza_labeled$CellLine.Treatment
mixed_Aza_labeled$CellLine.Treatment.Culture <- paste(Idents(mixed_Aza_labeled), mixed_Aza_labeled$Culture, sep = "_")
# Save correct Idents mixed_Aza_labeled
save(mixed_Aza_labeled, file = "/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/mixed_Aza_labeled_correct_Idents.RData")
save(mixed_Aza_labeled, file = "/Users/serena/scrnaseq/results/analysis/pureAzavcoAza/mixed_Aza_labeled_correct_Idents.RData")

# Integrate datasets: SNU719 and mixed_Aza--------------------------------------
mixed.SNU.Aza.list <- list(SNU719.Aza.april2017.seuratobject,SNU719.Aza.sept2017.seuratobject, mixed_Aza_labeled)
mixed.SNU.Aza.features <- SelectIntegrationFeatures(object.list = mixed.SNU.Aza.list, nfeatures = 3000)
mixed.SNU.Aza.list <- PrepSCTIntegration(object.list = mixed.SNU.Aza.list, anchor.features = mixed.SNU.Aza.features)
mixed.SNU.Aza.anchors <- FindIntegrationAnchors(object.list = mixed.SNU.Aza.list, normalization.method = "SCT", anchor.features = mixed.SNU.Aza.features)
mixed.SNU.Aza.integrated <- IntegrateData(anchorset = mixed.SNU.Aza.anchors, normalization.method = "SCT")
DefaultAssay(mixed.SNU.Aza.integrated) <- "integrated"

# save 
save(mixed.SNU.Aza.integrated, file = "/Users/serena/scrnaseq/results/analysis/pureAzavcoAza/mixed_SNU_Aza_integrated.Rdata")

# RunPCA Elbow Plot UMAP
mixed.SNU.Aza.integrated <- RunPCA(mixed.SNU.Aza.integrated, verbose = FALSE)
ElbowPlot(mixed.SNU.Aza.integrated)
PCAPlot(mixed.SNU.Aza.integrated,group.by = "sample")
mixed.SNU.Aza.integrated <- RunUMAP(mixed.SNU.Aza.integrated, reduction = "pca", dims = 1:15)

DimPlot(mixed.SNU.Aza.integrated, reduction = "umap") + ggtitle("Integrated Mono and Co SNU719 Aza Datasets")

DimPlot(mixed.SNU.Aza.integrated, group.by = "CellLine.Culture",
        cols = c('SNU719_Co' = 'tomato2', 'SNU719_Mono' = 'lightpink1',
                 'NCC24_Co' = "gray70",'1_Co' = "gray70",'2_Co' = "gray70",'4_Co' = "gray70",'5_Co' = "gray70",'6_Co' = "gray70",'7_Co' = "gray70")) +
  ggtitle("Integrated Mono and Co SNU719 Aza Datasets")

DimPlot(mixed.SNU.Aza.integrated, group.by = "CellLine",
        cols = c('SNU719' = 'tomato2',
                 'NCC24' = "forestgreen",'1' = "gray70",'2' = "gray70",'4' = "gray70",'5' = "gray70",'6' = "gray70",'7' = "gray70")) +
  ggtitle("Integrated Mono and Co SNU719 Aza Datasets")

Idents(mixed.SNU.Aza.integrated) <- "CellLine"
DimPlot(mixed.SNU.Aza.integrated, split.by = "Culture")


# Save integrated seurat object
save(mixed.SNU.Aza.integrated, file="/Users/serena/scrnaseq/results/analysis/pureAzavcoAza/integrate_mixed_SNU_Aza_withUMAP.RData")

# Set-up SNU object for DE: Normalize, Set Idents, ScaleData  ----------------------------------------------------------
# Set assay to RNA
DefaultAssay(mixed.SNU.Aza.integrated) <- "RNA"
# Normalize RNA data for visualization purposes
mixed.SNU.Aza.integrated <- NormalizeData(mixed.SNU.Aza.integrated, verbose = FALSE)
# Scale Data, and save the object again
# FVF to plot 
mixed.SNU.Aza.integrated <- FindVariableFeatures(mixed.SNU.Aza.integrated, selection.method = "vst", nfeatures = "3000")
# top10 <- head(VariableFeatures(mixed.SNU.Aza.integrated), 10)
# p3 <- VariableFeaturePlot(mixed.SNU.Aza.integrated)
# LabelPoints(plot = p3, points = top10, repel = TRUE)
mixed.SNU.Aza.integrated <- ScaleData(mixed.SNU.Aza.integrated, 
                                       features = VariableFeatures(mixed.SNU.Aza.integrated), 
                                       vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
# Save normalized, scaled object
save(mixed.SNU.Aza.integrated, file="/Users/serena/scrnaseq/results/analysis/pureAzavcoAza/integrate_mixed_SNU_Aza_withUMAP_and_ScaleData.RData")

# FCM: SNU Across Culture--------------------------------------
DefaultAssay(SNU719.Culture.Cells) <- "RNA"
Idents(SNU719.Culture.Cells) <- "CellLine"
# #Find conserved genes for SNU719-Aza across culture (Co v Mono)
SNU719.markers <- FindConservedMarkers(SNU719.Culture.Cells, ident.1 = 'SNU719', grouping.var = "Culture", verbose = FALSE)
head(SNU719.markers)

# Explore the marker genes for  SNU719 in FP
# Marker genes: High avglog_FC, pct1 >> pct2.
FeaturePlot(mixed.SNU.Aza.integrated, features = c("SPINK4"))

Idents(mixed.SNU.Aza.integrated) <- "CellLine"
markers.to.plot <- c("HIST1H4C","CCK","MT2A","LYZ")
DotPlot(mixed.SNU.Aza.integrated, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "Culture") + RotatedAxis()


# Average Expression of SNU719-----------------------------------------------------------
Idents(mixed.SNU.Aza.integrated) <- mixed.SNU.Aza.integrated$CellLine
SNU719.cells <- subset(mixed.SNU.Aza.integrated, idents = c("SNU719"))
Idents(SNU719.cells) <-mixed.SNU.Aza.integrated$Culture
avg.SNU719.cells <- log1p(AverageExpression(SNU719.cells, verbose = FALSE)$RNA)
p1 <- ggplot(avg.SNU719.cells, aes(`Mono`, `Co`)) + geom_point() + ggtitle("Aza-treated SNU719 in Mono vs Aza-treated SNU719 in Co-Culture")
HoverLocator(plot = `p1`)
genes.to.label = c("IFI6", "IFI27", "STAT1", "NOP53", "SLC25A6", "PRDX11","ACB","TPD52L1")
LabelPoints(p1, points = genes.to.label)
#DE: SNU--------------------------------------
# Edit Idents 
# mixed.SNU.Aza.integrated$CellLine.Treatment.Culture <- paste(mixed.SNU.Aza.integrated$CellLine.Treatment, mixed.SNU.Aza.integrated$Culture, sep = "_" )

# Load annotations
annotations <- read.csv("/users/serena/scrnaseq/data/annotation.csv")

# Create new object to perform analysis on:
mixed.SNU.Aza.integrated.new <- mixed.SNU.Aza.integrated
Idents(mixed.SNU.Aza.integrated.new)<-mixed.SNU.Aza.integrated.new$CellLine.Culture
DefaultAssay(mixed.SNU.Aza.integrated.new) <- "RNA"

# FindMarkers
Aza.SNU719.Culture.response <- FindMarkers(mixed.SNU.Aza.integrated.new, ident.1 = "SNU719_Co", ident.2 = "SNU719_Mono", verbose = TRUE, min.diff.pct = 0.3)
Aza.SNU719.Culture.response <- Aza.SNU719.Culture.response %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

# save list of de genes
write.csv(Aza.SNU719.Culture.response,file = "/Users/serena/scrnaseq/results/analysis/pureAzavcoAza/All_DE_genes_SNU719_pureAzavcoAza.csv")
save(Aza.SNU719.Culture.response, file = "/Users/serena/scrnaseq/results/analysis/pureAzavcoAza/All_DE_genes_SNU719_pureAzavcoAza.Rdata")

# Visualizing DE Genes  -----------------------------------------------------------
# List top 20 up reg and top 20 down reg
top20.Aza.SNU719.Culture.markers <- (top_n(Aza.SNU719.Culture.response, n = 20, wt = avg_logFC))$gene
top20.Aza.SNU719.Culture.markers <- append(top20.Aza.SNU719.Culture.markers, (top_n(Aza.SNU719.Culture.response, n = 20, wt = -avg_logFC))$gene)

# Create object of just SNU719_Mono and SNU719_Co cells
Idents(mixed.SNU.Aza.integrated.new)<-mixed.SNU.Aza.integrated.new$CellLine.Culture

Aza.SNU719.Cells <- subset(mixed.SNU.Aza.integrated.new, idents = c("SNU719_Mono", "SNU719_Co"))

# Plot Heatmap of these genes
DefaultAssay(Aza.SNU719.Cells) <- "RNA"
DoHeatmap(Aza.SNU719.Cells, group.by = "CellLine.Treatment.Culture",features = top20.Aza.SNU719.Culture.markers, 
          group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
  scale_fill_distiller(palette = "RdBu")

# Feature Plot of selected genes across treatment condition
# modify genes under 'features'
# moddify object to be plotted as desiredzw
FeaturePlot(Aza.SNU719.Cells, features = c("CDKN2A","ID2","IGFBP3"), 
            split.by = "Treatment", max.cutoff = 3, cols = c("grey", "red"))

# ,"TP53", "CDH1", "PTEN"

# Violin plot of selected genes across treatment condition
# modify genes under 'features'
plots <- VlnPlot(Aza.SNU719.Cells, 
                 features = c("IFI27","SKP1","TSPAN8" ,"ID1"), 
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




# Break --------------------------
# Integrate and normalize + scale datasets: SNU x2, mixed_Aza  -----------------------------------------------------------
mixed.SNU.Aza.list <- list(SNU719.Aza.april2017.seuratobject,SNU719.Aza.sept2017.seuratobject, mixed_Aza_labeled)
mixed.SNU.Aza.features <- SelectIntegrationFeatures(object.list = mixed.SNU.Aza.list, nfeatures = 3000)
mixed.SNU.Aza.list <- PrepSCTIntegration(object.list = mixed.SNU.Aza.list, anchor.features = mixed.SNU.Aza.features)
mixed.SNU.Aza.anchors <- FindIntegrationAnchors(object.list = mixed.SNU.Aza.list, normalization.method = "SCT", anchor.features = mixed.SNU.Aza.features)
mixed.SNU.Aza.integrated <- IntegrateData(anchorset = mixed.SNU.Aza.anchors, normalization.method = "SCT")

# Set assay
DefaultAssay(mixed.SNU.Aza.integrated) <- "integrated"

# save 
save(mixed.SNU.Aza.integrated, file = "mixed_SNU_Aza_integrated.Rdata")

# RunPCA Elbow Plot UMAP
mixed.SNU.Aza.integrated <- RunPCA(mixed.SNU.Aza.integrated, verbose = FALSE)
png("Elbow Plot Integrating mixedAza and pureSNU-Aza datasets.png")
ElbowPlot(mixed.SNU.Aza.integrated)
dev.off()
png("PCA Plot Integrating mixedAza and pureSNU-Aza datasets.png")
PCAPlot(mixed.SNU.Aza.integrated,group.by = "sample")
dev.off()
mixed.SNU.Aza.integrated <- RunUMAP(mixed.SNU.Aza.integrated, reduction = "pca", dims = 1:15)

table(Idents(mixed.SNU.Aza.integrated))

DimPlot(mixed.SNU.Aza.integrated, reduction = "umap") + ggtitle("Integrated Mono and Co SNU719 Aza Datasets")

DimPlot(mixed.SNU.Aza.integrated, group.by = "CellLine.Culture",
        cols = c('SNU719_Co' = 'tomato2', 'SNU719_Mono' = 'lightpink1',
                 'NCC24_Co' = "gray70",'1_Co' = "gray70",'2_Co' = "gray70",'4_Co' = "gray70",'5_Co' = "gray70",'6_Co' = "gray70",'7_Co' = "gray70")) +
  ggtitle("Integrated Mono and Co SNU719 Aza Datasets")
DimPlot(mixed.SNU.Aza.integrated, group.by = "CellLine",
        cols = c('SNU719' = 'tomato2',
                 'NCC24' = "forestgreen",'1' = "gray70",'2' = "gray70",'4' = "gray70",'5' = "gray70",'6' = "gray70",'7' = "gray70")) +
  ggtitle("Integrated Mono and Co SNU719 Aza Datasets")

Idents(mixed.SNU.Aza.integrated) <- "CellLine"
DimPlot(mixed.SNU.Aza.integrated, split.by = "Culture")

# Save integrated seurat object
save(mixed.SNU.Aza.integrated, file="/Users/serena/scrnaseq/results/analysis/pureAzavcoAza/integrate_mixed_SNU_Aza_withUMAP.RData")

# Set-up SNU object for DE: Normalize, Set Idents, ScaleData  ----------------------------------------------------------
# Set assay to RNA
DefaultAssay(mixed.SNU.Aza.integrated) <- "RNA"
# Normalize RNA data for visualization purposes
mixed.SNU.Aza.integrated <- NormalizeData(mixed.SNU.Aza.integrated, verbose = FALSE)
mixed.SNU.Aza.integrated <- FindVariableFeatures(mixed.SNU.Aza.integrated, selection.method = "vst", nfeatures = 3000)
# Scale Data, and save the object again
mixed.SNU.Aza.integrated <- ScaleData(mixed.SNU.Aza.integrated, vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))

# Save normalized, scaled object
save(mixed.SNU.Aza.integrated, file="/Users/serena/scrnaseq/results/analysis/pureAzavcoAza/integrate_mixed_SNU_Aza_withUMAP_and_ScaleData.RData")

#DE: SNU--------------------------------------
# Load annotations
annotations <- read.csv("/users/serena/scrnaseq/data/annotation.csv")
# Create new object to perform analysis on:
mixed.SNU.Aza.integrated.new <- mixed.SNU.Aza.integrated
# Set Idents - What idents to set to?
# I am testing SNU719_Co vs SNU719_Mono
Idents(mixed.SNU.Aza.integrated.new)<-mixed.SNU.Aza.integrated.new$CellLine.Culture
# Set RNA assay
DefaultAssay(mixed.SNU.Aza.integrated.new) <- "RNA"

# FindMarkers
SNU719.Culture.response <- FindMarkers(mixed.SNU.Aza.integrated.new, ident.1 = "SNU719_Co", ident.2 = "SNU719_Mono", verbose = TRUE, min.diff.pct = 0.3)
SNU719.Culture.response <- SNU719.Culture.response %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

# List top 20 up reg and top 20 down reg
genes.SNU719.Culture <- (top_n(SNU719.Culture.response, n = 20, wt = avg_logFC))$gene
genes.SNU719.Culture <- append(genes.SNU719.Culture, (top_n(SNU719.Culture.response, n = 20, wt = -avg_logFC))$gene)

# Create object of just SNU719_Aza and SNU719_Culture cells
SNU719.Culture.Cells <- subset(mixed.SNU.Aza.integrated.new, idents = c("SNU719_Mono", "SNU719_Co"))

# Plot Heatmap of these genes
DefaultAssay(SNU719.Culture.Cells) <- "RNA"
DoHeatmap(SNU719.Culture.Cells, group.by = "CellLine.Culture",features = genes.SNU719.Culture, group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu")
