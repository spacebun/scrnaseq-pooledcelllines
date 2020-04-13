#!/usr/bin/env Rscript

# README ------------------------------------------------------------------

# How does being in co-culture affect the treatment response (gene expression changes) of SNU719 cells?

# The lists of 40 (top 20 up- and top 20 down-regulated) differentially expressed (DE) genes for
# mono-cultured SNU719 in DMSO vs Aza and
# co-cultured SNU719 in DMSO vs Aza 
# are examined for overlapping genes using a Venn Diagram (https://bioinfogp.cnb.csic.es/tools/venny/)

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
setwd("/Users/serena/scrnaseq/results/analysis/comparetreatmentacrossculture")

# Load the integrated pure SNU719 (DMSO vs Aza)
load(file = "/Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza/integrate_mixed_scaledRNAassay.RData")
load(file = "/Users/serena/scrnaseq/results/analysis/mixedDMSOvmixedAza/All_DE_genes_SNU719_mixedDMSOvmixedAza.Rdata")

# Load the integrated mixed SNU719 (DMSO vs Aza)
load(file = "/Users/serena/scrnaseq/results/analysis/pureDMSOvpureAza/integrate_SNU_scaledRNAassay.RData")
load(file = "/Users/serena/scrnaseq/results/analysis/pureDMSOvpureAza/All_DE_genes_SNU719_pureDMSOvpureAza.Rdata")


# Create object of just SNU719_DMSO and SNU719_Aza cells from co-cultured dataset
Idents(mixed.integrated) <-mixed.integrated$CellLine.Treatment
mixed.SNU719 <- subset(mixed.integrated, idents = c("SNU719_DMSO", "SNU719_Aza"))

# subset top 20 DE genes for each treatment condition
# mono
top20.pure.SNU719.Treatment.markers <- pure.SNU719.Treatment.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
# co 
top20.mixed.SNU719.Treatment.markers  <- (top_n(mixed.SNU719.Treatment.markers, n = 20, wt = avg_logFC))$gene
top20.mixed.SNU719.Treatment.markers <- append(top20.mixed.SNU719.Treatment.markers, (top_n(mixed.SNU719.Treatment.markers, n = 20, wt = -avg_logFC))$gene)

# Using Venny, obtain list of genes categorzed by which dataset they are exclusive to
de_genes_common_list<-read.csv(file  = "de_genes_common_list.csv")
list <- as.character(de_genes_common_list$gene)

# Visualize results on individually normalized objects-----------------------------------------------------------
# plot list of all genes on mono 
DoHeatmap(SNU.integrated, group.by = "CellLine.Treatment.Culture",features = list, group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
  scale_fill_distiller(palette = "RdBu")
# plot list of all genes on co 
DoHeatmap(mixed.SNU719, group.by = "CellLine.Treatment.Culture",features = list, group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
  scale_fill_distiller(palette = "RdBu")

# Common genes
common.genes <- c("PHGR1", "PTMA",
"AC020656.1",
"HIST1H4C",
"ADIRF",
"ANXA2",
"ARPC1B",
"LGALS1",
"CCK",
"CEACAM6")
# plot common genes on mono 
DoHeatmap(SNU.integrated, group.by = "CellLine.Treatment.Culture",features = common.genes, group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu")
# plot common genes on co
DoHeatmap(mixed.SNU719, group.by = "CellLine.Treatment.Culture",features = common.genes, group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 1, angle = 0) +
  scale_fill_distiller(palette = "RdBu")
# Merge all datasets and normalize again-----------------------------------------------------------
load("/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/mixed_DMSO_labeled_v0.RData")
load("/Users/serena/scrnaseq/results/unmixing-v1/use-merge-as-ref/v0/mixed_Aza_labeled_v0.RData")
load("/Users/serena/scrnaseq/results/pure_indiv/april2017/SNU719-Aza-april2017/v1/SCT-SNU719-Aza-april2017-v1.RData")
load("/Users/serena/scrnaseq/results/pure_indiv/april2017/SNU719-DMSO-april2017/v1/SCT-SNU719-DMSO-april2017-v1.RData")
load("/Users/serena/scrnaseq/results/pure_indiv/july2017/SNU719-DMSO-july2017/v1/SCT-SNU719-DMSO-july2017-v1.RData")
load("/Users/serena/scrnaseq/results/pure_indiv/sept2017/SNU719-Aza-sept2017/v1/SCT-SNU719-Aza-sept2017-v1.RData")

# merge all datasets
ctac.merged<- merge(x = SNU.integrated,
                         y = c(mixed.integrated), 
                    add.cell.ids = c("pure.SNU.integrated", "mixed.SNU.integrated"),
                    project = "ctac_merged")
Idents(ctac.merged) <- ctac.merged$CellLine.Treatment.Culture
DefaultAssay(ctac.merged) <- "RNA"
ctac.merged <- subset(ctac.merged, idents = c("SNU719_DMSO_Co", "SNU719_DMSO_Mono", "SNU719_Aza_Co", "SNU719_Aza_Mono"))
ctac.merged <- NormalizeData(ctac.merged)
ctac.merged <- FindVariableFeatures(ctac.merged, nfeatures = 3000)
all.genes <- rownames(ctac.merged)
ctac.merged <- ScaleData(ctac.merged, features = all.genes, vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
save(ctac.merged, file = "/Users/serena/scrnaseq/results/analysis/comparetreatmentacrossculture/ctac_merged.Rdata")
DoHeatmap(ctac.merged, group.by = "CellLine.Treatment.Culture",features = list, group.bar = TRUE, draw.lines = TRUE, size = 4, hjust = 0.5, angle = 0) +
  scale_fill_distiller(palette = "RdBu")
