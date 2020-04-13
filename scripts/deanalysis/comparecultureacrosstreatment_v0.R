#!/usr/bin/env Rscript

# README ------------------------------------------------------------------
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

# Set working directory
setwd("/Users/serena/scrnaseq/results/analysis/comparecultureacrosstreatment")

# Load the integrated DMSO SNU719 (Co vs Mono)
load(file = "/Users/serena/scrnaseq/results/analysis/pureDMSOvcoDMSO/integrate_mixed_SNU_DMSO_withUMAP_and_ScaleData.RData")
load(file = "/Users/serena/scrnaseq/results/analysis/pureDMSOvcoDMSO/All_DE_genes_SNU719_pureDMSOvcoDMSO.Rdata")

# Load the integrated Aza SNU719 (Co vs Mono)
load(file = "/Users/serena/scrnaseq/results/analysis/pureAzavcoAza/integrate_mixed_SNU_Aza_withUMAP_and_ScaleData.RData")
load(file = "/Users/serena/scrnaseq/results/analysis/pureAzavcoAza/All_DE_genes_SNU719_pureAzavcoAza.Rdata")

# subset top 20 DE genes for each culture condition
# DMSO
top20.DMSO.SNU719.Culture.markers <- (top_n(DMSO.SNU719.Culture.response, n = 20, wt = avg_logFC))$gene
top20.DMSO.SNU719.Culture.markers <- append(top20.DMSO.SNU719.Culture.markers, (top_n(DMSO.SNU719.Culture.response, n = 20, wt = -avg_logFC))$gene)
write.csv(top20.DMSO.SNU719.Culture.markers, file = "top 20 markers DMSO across culture.csv")

# Aza 
top20.Aza.SNU719.Culture.markers <- (top_n(Aza.SNU719.Culture.response, n = 20, wt = avg_logFC))$gene
top20.Aza.SNU719.Culture.markers <- append(top20.Aza.SNU719.Culture.markers, (top_n(Aza.SNU719.Culture.response, n = 20, wt = -avg_logFC))$gene)
write.csv(top20.Aza.SNU719.Culture.markers, file = "top 20 markers Aza across culture.csv")

# Using Venny, obtain list of genes categorzed by which dataset they are exclusive to
de_genes_common_list<-read.csv(file  = "de_genes_common_list.csv")
list <- as.character(de_genes_common_list$gene)

# Visualize results on individually normalized objects-----------------------------------------------------------
# plot list of all genes on DMSO 
Idents(mixed.SNU.DMSO.integrated) <- mixed.SNU.DMSO.integrated$CellLine.Culture
DefaultAssay(mixed.SNU.DMSO.integrated) <- "RNA"
DMSO.SNU719.Cells <- subset(mixed.SNU.DMSO.integrated, idents = c("SNU719_Mono", "SNU719_Co"))
DoHeatmap(DMSO.SNU719.Cells, group.by = "CellLine.Treatment.Culture",features = list, group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
  scale_fill_distiller(palette = "RdBu")
# plot list of all genes on Aza 
Idents(mixed.SNU.Aza.integrated) <- mixed.SNU.Aza.integrated$CellLine.Culture
DefaultAssay(mixed.SNU.Aza.integrated) <- "RNA"
Aza.SNU719.Cells <- subset(mixed.SNU.Aza.integrated, idents = c("SNU719_Mono", "SNU719_Co"))
DoHeatmap(Aza.SNU719.Cells, group.by = "CellLine.Treatment.Culture",features = list, group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
  scale_fill_distiller(palette = "RdBu")

# Merge all datasets and normalize again-----------------------------------------------------------
ccat.merged<- merge(x = mixed.SNU.DMSO.integrated,
                    y = c(mixed.SNU.Aza.integrated),
                    project = "ccat_merged", 
                    add.cell.id = c("mixed.SNU.DMSO.integrated", "mixed.SNU.Aza.integrated"))
Idents(ccat.merged) <- ccat.merged$CellLine.Treatment.Culture
DefaultAssay(ccat.merged) <- "RNA"
ccat.merged <- subset(ccat.merged, idents = c("SNU719_DMSO_Mono", "SNU719_DMSO_Co", "SNU719_Aza_Mono", "SNU719_Aza_Co"))
ccat.merged <- NormalizeData(ccat.merged)
ccat.merged <- FindVariableFeatures(ccat.merged, nfeatures = 3000)
all.genes <- rownames(ccat.merged)
ccat.merged <- ScaleData(ccat.merged, features = all.genes, vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
save(ccat.merged, file = "/Users/serena/scrnaseq/results/analysis/comparecultureacrosstreatment/ccat.merged.Rdata")

DoHeatmap(ccat.merged, group.by = "CellLine.Treatment.Culture",features = list, group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
  scale_fill_distiller(palette = "RdBu") 

plots <- VlnPlot(ccat.merged, 
                 features = c("IFI6", "RPS17", "IFI27", "MGST2", "SKP1", "CCK"), 
                 split.by = "CellLine.Culture", 
                 group.by = 'CellLine.Treatment',
                 pt.size = 0, 
                 combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x +
                  theme(axis.title = element_text(size=8),
                        axis.text.x = element_text(size=8,angle = 0, hjust = 0.5),
                        axis.text.y = element_text(size=8),
                        plot.title=element_text(size=8,face="bold"),
                        legend.position = 'top'))
CombinePlots(plots = plots, nrow = 3)

ccat.merged <- FindVariableFeatures(ccat.merged, nfeatures = 3000)
ccat.merged <- RunPCA(ccat.merged, features = VariableFeatures(ccat.merged))
ccat.merged <- RunUMAP(ccat.merged, reduction = 'pca', dims = 1:10)
DimPlot(ccat.merged, reduction = "umap", group.by = "CellLine.Treatment.Culture")
FeaturePlot(ccat.merged, features = c("TFF1"), 
            split.by = "Culture", 
            max.cutoff = 3, cols = c("grey", "red"))


# # merge all datasets
# ccat.merged<- merge(x = mixed.SNU.DMSO.integrated,
#                     y = c(mixed.SNU.Aza.integrated),
#                     project = "ccat_merged", 
#                     add.cell.id = c("mixed.SNU.DMSO.integrated", "mixed.SNU.Aza.integrated"))
# Idents(ccat.merged) <- ccat.merged$CellLine.Treatment.Culture
# DefaultAssay(ccat.merged) <- "RNA"
# ccat.merged <- NormalizeData(ccat.merged)
# ccat.merged <- FindVariableFeatures(ccat.merged, nfeatures = 3000)
# ccat.merged <- ScaleData(ccat.merged, features = VariableFeatures(ccat.merged), vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
# save(ccat.merged, file = "/Users/serena/scrnaseq/results/analysis/comparecultureacrosstreatment/ccat.merged.Rdata")
# SNU.ccat.merged <- subset(ccat.merged, idents = c("SNU719_DMSO_Mono", "SNU719_DMSO_Co", "SNU719_Aza_Mono", "SNU719_Aza_Co"))
# SNU.ccat.merged <- NormalizeData(SNU.ccat.merged)
# SNU.ccat.merged <- FindVariableFeatures(SNU.ccat.merged, nfeatures = 3000)
# SNU.ccat.merged <- ScaleData(SNU.ccat.merged, features = VariableFeatures(SNU.ccat.merged), vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
# DoHeatmap(ccat.merged, group.by = "CellLine.Treatment.Culture",features = list, group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
#   scale_fill_distiller(palette = "RdBu") 
