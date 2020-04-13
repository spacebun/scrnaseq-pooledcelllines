#!/usr/bin/env Rscript

# README ------------------------------------------------------------------
# Sample Name: mixed_DMSO
# Goal: Create and pre-process sample as Seurat object,
# to be integrated with other samples

# Script version: 1
# Updates: Replace filtering step: Use SCTransform instead of subset to filter out cells.
# Removed quality metric plots. Edited steps in SCtransform in regressing out cell cycle.
# Newly saved SCtransformed object has RunPCA applied to it.

# Set-up packages, working directory, dataset -----------------------------------------------------------
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
library(beepr)

## Set directory
setwd("/Users/serena/scrnaseq/results/mixed_indiv/mixed_DMSO/v0")

## Load the Dataset
mixed.DMSO.data <- Read10X(data.dir = "/users/serena/scrnaseq/results/cellranger-outs-mixed/mixed_DMSO/filtered_feature_bc_matrix")

## Initialize Seurat object with the raw (non-normalized Data).
# note project name fills in orig.ident slot in object
mixed.DMSO.seuratobject <- CreateSeuratObject(counts = mixed.DMSO.data, project = "mixed.DMSO", min.cells = 3, min.features = 500)

# Generating quality metrics and selecting (filtering out) cells  -----------------------------
# Record Information on Seurat object, dense.size, sparse.size and dense.size/sparse.size
cat("mixed_DMSO-Seurat Object", file = "mixed_DMSO-seuratresults.txt")
cat("\n", file = "mixed_DMSO-seuratresults.txt", append = TRUE)
capture.output(mixed.DMSO.seuratobject,file = "mixed_DMSO-seuratresults.txt", append = TRUE)

dense.size <- object.size(as.matrix(mixed.DMSO.data))
cat("\n\n", file = "mixed_DMSO-seuratresults.txt", append = TRUE)
cat("dense.size", file = "mixed_DMSO-seuratresults.txt", append = TRUE)
cat("\n", file = "mixed_DMSO-seuratresults.txt", append = TRUE)
capture.output(dense.size,file = "mixed_DMSO-seuratresults.txt", append = TRUE)

sparse.size <- object.size(mixed.DMSO.data)
cat("\n\n", file = "mixed_DMSO-seuratresults.txt", append = TRUE)
cat("sparse.size", file = "mixed_DMSO-seuratresults.txt", append = TRUE)
cat("\n", file = "mixed_DMSO-seuratresults.txt", append = TRUE)
capture.output(sparse.size,file = "mixed_DMSO-seuratresults.txt", append = TRUE)

cat("\n\n", file = "mixed_DMSO-seuratresults.txt", append = TRUE)
cat("dense.size/sparse.size", file = "mixed_DMSO-seuratresults.txt", append = TRUE)
cat("\n", file = "mixed_DMSO-seuratresults.txt", append = TRUE)
capture.output(dense.size/sparse.size,file = "mixed_DMSO-seuratresults.txt", append = TRUE)

# Add number of genes per UMI for each cell to metadata
mixed.DMSO.seuratobject$log10GenesPerUMI <- log10(mixed.DMSO.seuratobject$nFeature_RNA) / log10(mixed.DMSO.seuratobject$nCount_RNA)

## Calculate percentage of counts originating from the set of features with MT- as set of mitochondrial genes
# exclude cells that are damaged by checking relative expression of mitochondrially derived genes
mixed.DMSO.seuratobject[["percent.mt"]] <- PercentageFeatureSet(mixed.DMSO.seuratobject, pattern = "^MT-")

# Create metadata dataframe as variable to insert addtl metrics for QC without affecting OG merged object
metadata <- mixed.DMSO.seuratobject@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Add sample names
metadata$sample <- NA
metadata$sample<- "mixed_DMSO"

# Add metadata back to Seurat object
mixed.DMSO.seuratobject@meta.data <- metadata

# Add more information to dataset
mixed.DMSO.seuratobject$CellLine <- 'Mixed'
mixed.DMSO.seuratobject$Status <- 'Cancer'
mixed.DMSO.seuratobject$Subtype <- 'Mixed'
mixed.DMSO.seuratobject$Treatment <- 'DMSO'
mixed.DMSO.seuratobject$Date <- 'NA'
mixed.DMSO.seuratobject$TP53Mutation <- 'Mixed'
mixed.DMSO.seuratobject$Gender <- 'Mixed'
mixed.DMSO.seuratobject$Culture <- 'Co'

## Assess Quality Metric Plots (save each plot)
dir.create("Quality Metric Plots")
setwd("/Users/serena/scrnaseq/results/mixed_indiv/mixed_DMSO/v0/Quality Metric Plots")

# Visualize UMI counts (transcripts) per cell
# Visualize the number UMIs/transcripts per cell
png('nUMIs per cell.png')
metadata %>%
   ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) +
   geom_density(alpha = 0.2) +
   scale_x_log10() +
   theme_classic() +
   ylab("Cell density") +
   geom_vline(xintercept = 10000)+
   ggtitle("nUMIs per cell")
dev.off()

# Visualize the distribution of genes detected per cell via histogram
png('Distn Genes per cell.png')
metadata %>%
 ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) +
 geom_density(alpha = 0.2) +
 theme_classic() +
 scale_x_log10() +
 #geom_vline(xintercept = 2800)+
 ggtitle("Genes per cell")
dev.off()

# Visualize the distribution of genes detected per cell via boxplot
png('Distn log10 genes per cell boxplot.png')
metadata %>%
  ggplot(aes(x=sample, y=log10(nFeature_RNA), fill=sample)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Genes per cell")
dev.off()

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
png('nGenes vs nUMIs.png')
metadata %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
  geom_point() +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  # geom_vline(xintercept = 9000) +
  # geom_hline(yintercept = 3100) +
  facet_wrap(~sample) +
  ggtitle("nGenes vs nUMIs")
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per cell
png('Distn Mito Genes per cell.png')
metadata %>%
   ggplot(aes(color=sample, x=percent.mt, fill=sample)) +
   geom_density(alpha = 0.2) +
   scale_x_log10() +
   theme_classic() +
   geom_vline(xintercept = 1.0) +
   ggtitle("Mitochondrial Genes per Cell")
dev.off()

# Visualize the overall novelty of the gene expression by visualizing the genes detected per UMI
png('Novelty plot.png')
metadata %>%
 ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
 geom_density(alpha = 0.2) +
 theme_classic() +
 geom_vline(xintercept = 0.875) +
 ggtitle("Genes detected per UMI")
dev.off()
 # VlnPlot(mixed.DMSO.seuratobject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png('VlnPlot Genes per cell.png')
VlnPlot(mixed.DMSO.seuratobject, features = c("nFeature_RNA") ) + ylab("Expression Level")
dev.off()
png('VlnPlot nCount_RNA per cell.png')
VlnPlot(mixed.DMSO.seuratobject, features = c("nCount_RNA") ) + ylab("Expression Level")
dev.off()
png('VlnPlot percentmt per cell.png')
VlnPlot(mixed.DMSO.seuratobject, features = c("percent.mt") ) + ylab("Expression Level") + scale_y_continuous(breaks=seq(0,100,5))
dev.off()
VlnPlot(mixed.DMSO.seuratobject, features = c("percent.mt"),  y.max = 100) + 
  ggtitle("Distribution of cells according to the percent of metochondrial reads in cell") +
  ylab("Expression Level") + 
  scale_y_continuous(breaks=seq(0,100,5)) + 
  theme(axis.title = element_text(size=10),axis.text.x = element_text(size=10,angle = 0, hjust = 0.5),axis.text.y = element_text(size=10),plot.title=element_text(size=10,face="bold")) +   
  # geom_hline(yintercept = 20) + 
  theme(legend.position = 'none')

# visualize feature-feature relationships with FeatureScatter
# Save manually
plot1 <- FeatureScatter(mixed.DMSO.seuratobject, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mixed.DMSO.seuratobject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# # Filter out low quality reads using selected thresholds----------------------
mixed.DMSO.seuratobject <- subset(mixed.DMSO.seuratobject,
  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 &
    nCount_RNA < 10000 &
   percent.mt < 20)

# Run SCTransform: ----------------------
# Normalizes data, Scales Data, Finds variable features
# and regress out percent.mt etc
# SCTransform operate on original unnormalized data in 'RNA' assay slot
#  and saves normalized data in 'SCT' slot
mixed.DMSO.seuratobject <-SCTransform(mixed.DMSO.seuratobject,
  assay = 'RNA',
  new.assay.name = 'SCT',
  vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'), verbose = FALSE)

# Saving filtered object
setwd("/Users/serena/scrnaseq/results/mixed_indiv/mixed_DMSO/v0")
save(mixed.DMSO.seuratobject, file="filtered-mixed_DMSO-v0.RData")

# Score cells for cell cycle
# CellCycleScoring operates on SCT assay slot.
load("/users/serena/scrnaseq/data/cycle.rda")
mixed.DMSO.seuratobject <- CellCycleScoring(mixed.DMSO.seuratobject,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes,
                                 assay = 'SCT',
                                 set.ident = TRUE)

# Perform PCA to check cell cycle effect
mixed.DMSO.seuratobject <- RunPCA(mixed.DMSO.seuratobject)

# Plot the PCA colored by cell cycle phase
# Should observe variation by cell cycle
png('Cell Cycle PCA.png')
DimPlot(mixed.DMSO.seuratobject,
       reduction = "pca",
       group.by = "Phase")
dev.off()

# Run SCTransform:
# This time, operate on original unnormalized data in 'RNA' assay slot and
# save to 'SCT' slot, overwriting previous SCT.
mixed.DMSO.seuratobject <-SCTransform(mixed.DMSO.seuratobject,
  assay = 'RNA',
  new.assay.name = 'SCT',
  vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))

# Perform PCA
mixed.DMSO.seuratobject <- RunPCA(mixed.DMSO.seuratobject, features = VariableFeatures(object = mixed.DMSO.seuratobject))

# Now, a PCA on the variable genes no longer returns components associated with cell cycle
setwd("/Users/serena/scrnaseq/results/mixed_indiv/mixed_DMSO/v0")
png('Cell Cycle Regressed PCA.png')
DimPlot(mixed.DMSO.seuratobject,
        reduction = "pca",
        group.by= "Phase")
dev.off()

# Saving SCtransformed object---------------------------------------------------------
# To be used as input to integration. RunPCA applied,
setwd("/Users/serena/scrnaseq/results/mixed_indiv/mixed_DMSO/v0")
save(mixed.DMSO.seuratobject, file="SCT-mixed_DMSO-v0.RData")

# Linear Dimensional Reduction and checking Cell Cycle------------------------------------

# ## Examine and visualize PCA results a few different ways
# cat("\n\n", file = "mixed_DMSO-seuratresults.txt", append = TRUE)
# cat("Most variable genes driving PCs", file = "mixed_DMSO-seuratresults.txt", append = TRUE)
# cat("\n", file = "mixed_DMSO-seuratresults.txt", append = TRUE)
# capture.output(print(mixed.DMSO.seuratobject[["pca"]], dims = 1:10, nfeatures = 5),file = "mixed_DMSO-seuratresults.txt", append = TRUE)

pdf('mixed_DMSO-PCAPlots.pdf')
VizDimLoadings(mixed.DMSO.seuratobject, dims = 1:2, reduction = "pca") # all cells plotted on first two PC
DimPlot(mixed.DMSO.seuratobject, reduction = "pca", group.by = "orig.ident")
DimHeatmap(mixed.DMSO.seuratobject, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(mixed.DMSO.seuratobject, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

# Determine ‘dimensionality’ of the dataset  ----------------------
# The bend at the elbow plot at X PCs, suggests the majority of signal is captured in the first X PCs.
# Choose how many PC dimensions you want to include in the next step of clustering based on the elbow plot.
png('mixed_DMSO-ElbowPlot.png')
ElbowPlot(mixed.DMSO.seuratobject)
dev.off()
beep(2)
