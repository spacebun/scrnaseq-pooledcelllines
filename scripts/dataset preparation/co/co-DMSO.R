#!/usr/bin/env Rscript

# README ------------------------------------------------------------------
# Sample Name: co-DMSO
# Set-up packages, working directory, dataset -----------------------------------------------------------
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
library(beepr)

## Set directory
setwd("/scrnaseq/results/dataset preparation/co/co-DMSO")

## Load the Dataset
co.DMSO.data <- Read10X(data.dir = "/scrnaseq/results/cellranger-outs-mixed/mixed_DMSO/filtered_feature_bc_matrix")

## Initialize Seurat object with the raw (non-normalized Data).
# note project name fills in orig.ident slot in object
co.DMSO.seuratobject <- CreateSeuratObject(counts = co.DMSO.data, project = "co.DMSO", min.cells = 3, min.features = 500)

# Generating quality metrics and selecting (filtering out) cells  -----------------------------

# Add number of genes per UMI for each cell to metadata
co.DMSO.seuratobject$log10GenesPerUMI <- log10(co.DMSO.seuratobject$nFeature_RNA) / log10(co.DMSO.seuratobject$nCount_RNA)

## Calculate percentage of counts originating from the set of features with MT- as set of mitochondrial genes
# exclude cells that are damaged by checking relative expression of mitochondrially derived genes
co.DMSO.seuratobject[["percent.mt"]] <- PercentageFeatureSet(co.DMSO.seuratobject, pattern = "^MT-")

# Create metadata dataframe as variable to insert addtl metrics for QC without affecting OG merged object
metadata <- co.DMSO.seuratobject@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Add sample names
metadata$sample <- NA
metadata$sample<- "co_DMSO"

# Add metadata back to Seurat object
co.DMSO.seuratobject@meta.data <- metadata

# Add more information to dataset
co.DMSO.seuratobject$CellLine <- 'Co'
co.DMSO.seuratobject$Status <- 'Cancer'
co.DMSO.seuratobject$Subtype <- 'Co'
co.DMSO.seuratobject$Treatment <- 'DMSO'
co.DMSO.seuratobject$Date <- 'NA'
co.DMSO.seuratobject$TP53Mutation <- 'Co'
co.DMSO.seuratobject$Gender <- 'Co'
co.DMSO.seuratobject$Culture <- 'Co'
co.DMSO.seuratobject$CellLine.Treatment <- "Co_DMSO"
co.DMSO.seuratobject$CellLine.Culture <- "Co_Mono"
co.DMSO.seuratobject$CellLine.Treatment.Culture <- "Co_DMSO_Mono"

## Assess Quality Metric Plots (save each plot)
# Visualize UMI counts (transcripts) per cell

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
png('before QC nGenes vs nUMIs.png')
metadata %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
  geom_point() +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 9000) +
  geom_hline(yintercept = 3500) +
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
 # VlnPlot(co.DMSO.seuratobject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png('VlnPlot Genes per cell.png')
VlnPlot(co.DMSO.seuratobject, features = c("nFeature_RNA") ) + ylab("Expression Level")
dev.off()
png('VlnPlot nCount_RNA per cell.png')
VlnPlot(co.DMSO.seuratobject, features = c("nCount_RNA") ) + ylab("Expression Level")
dev.off()
png('VlnPlot percentmt per cell.png')
VlnPlot(co.DMSO.seuratobject, features = c("percent.mt") ) + ylab("Expression Level") + scale_y_continuous(breaks=seq(0,100,5))
dev.off()

png('Before QC VlnPlot percentmt per cell.png')
VlnPlot(co.DMSO.seuratobject, features = c("percent.mt"),  y.max = 100) +
  ggtitle("Distribution of cells according to the percent of metochondrial reads in cell") +
  ylab("Expression Level") +
  scale_y_continuous(breaks=seq(0,100,5)) +
  theme(axis.title = element_text(size=10),axis.text.x = element_text(size=10,angle = 0, hjust = 0.5),axis.text.y = element_text(size=10),plot.title=element_text(size=10,face="bold")) +
  geom_hline(yintercept = 20) +
  theme(legend.position = 'none')
dev.off()
# visualize feature-feature relationships with FeatureScatter
# Save manually
# plot1 <- FeatureScatter(co.DMSO.seuratobject, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(co.DMSO.seuratobject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))

# # Filter out low quality reads using selected thresholds----------------------
co.DMSO.seuratobject <- subset(co.DMSO.seuratobject,
  subset = nFeature_RNA > 200 & nFeature_RNA < 3500 &
    nCount_RNA < 9000 &
   percent.mt < 20)
metadata <- co.DMSO.seuratobject@meta.data

png('After QC nGenes vs nUMIs.png')
metadata %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
  geom_point() +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  # geom_vline(xintercept = 9000) +
  # geom_hline(yintercept = 3500) +
  facet_wrap(~sample) +
  ggtitle("nGenes vs nUMIs")
dev.off()
png('After QC VlnPlot percentmt per cell.png')
VlnPlot(co.DMSO.seuratobject, features = c("percent.mt"),  y.max = 100) +
  ggtitle("Distribution of cells according to the percent of metochondrial reads in cell") +
  ylab("Expression Level") +
  scale_y_continuous(limit=c(0,100)) +
  theme(axis.title = element_text(size=10),axis.text.x = element_text(size=10,angle = 0, hjust = 0.5),axis.text.y = element_text(size=10),plot.title=element_text(size=10,face="bold")) +
  geom_hline(yintercept = 20) +
  theme(legend.position = 'none')
dev.off()

# Run SCTransform: ----------------------
# Normalizes data, Scales Data, Finds variable features
# and regress out percent.mt etc
# SCTransform operate on original unnormalized data in 'RNA' assay slot
#  and saves normalized data in 'SCT' slot
co.DMSO.seuratobject <-SCTransform(co.DMSO.seuratobject,
  assay = 'RNA',
  new.assay.name = 'SCT',
  vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'), verbose = FALSE)


# Score cells for cell cycle
# CellCycleScoring operates on SCT assay slot.
load("/scrnaseq/data/cycle.rda")
co.DMSO.seuratobject <- CellCycleScoring(co.DMSO.seuratobject,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes,
                                 assay = 'SCT',
                                 set.ident = TRUE)

# Perform PCA to check cell cycle effect
co.DMSO.seuratobject <- RunPCA(co.DMSO.seuratobject)

# Plot the PCA colored by cell cycle phase
# Should observe variation by cell cycle
png('Cell Cycle PCA.png')
DimPlot(co.DMSO.seuratobject,
       reduction = "pca",
       group.by = "Phase")
dev.off()

# Run SCTransform:
# This time, operate on original unnormalized data in 'RNA' assay slot and
# save to 'SCT' slot, overwriting previous SCT.
co.DMSO.seuratobject <-SCTransform(co.DMSO.seuratobject,
  assay = 'RNA',
  new.assay.name = 'SCT',
  vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))

# Perform PCA
co.DMSO.seuratobject <- RunPCA(co.DMSO.seuratobject, features = VariableFeatures(object = co.DMSO.seuratobject))

# Now, a PCA on the variable genes no longer returns components associated with cell cycle
png('Cell Cycle Regressed PCA.png')
DimPlot(co.DMSO.seuratobject,
        reduction = "pca",
        group.by= "Phase")
dev.off()

# Saving SCtransformed object---------------------------------------------------------
# To be used as input to integration. RunPCA applied,
coculture.DMSO <- co.DMSO.seuratobject
save(coculture.DMSO, file="SCT-co_DMSO.RData")

# Linear Dimensional Reduction and checking Cell Cycle------------------------------------
pdf('co_DMSO-PCAPlots.pdf')
VizDimLoadings(co.DMSO.seuratobject, dims = 1:2, reduction = "pca") # all cells plotted on first two PC
DimPlot(co.DMSO.seuratobject, reduction = "pca", group.by = "orig.ident")
DimHeatmap(co.DMSO.seuratobject, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(co.DMSO.seuratobject, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

# Determine ‘dimensionality’ of the dataset  ----------------------
# The bend at the elbow plot at X PCs, suggests the majority of signal is captured in the first X PCs.
# Choose how many PC dimensions you want to include in the next step of clustering based on the elbow plot.
png('co_DMSO-ElbowPlot.png')
ElbowPlot(co.DMSO.seuratobject)
dev.off()
