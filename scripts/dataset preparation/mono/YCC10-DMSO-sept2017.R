#!/usr/bin/env Rscript

# README ------------------------------------------------------------------
# Sample Name: YCC10-DMSO-sept2017

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
setwd("/scrnaseq/results/dataset preparation/mono/sept2017/YCC10-DMSO-sept2017")

## Load the Dataset
YCC10.DMSO.sept2017.data <- Read10X(data.dir = "/scrnaseq/results/cellranger-outs/sept_2017/YCC10-DMSO-sept2017/filtered_feature_bc_matrix")

## Initialize Seurat object with the raw (non-normalized Data).
# note project name fills in orig.ident slot in object
YCC10.DMSO.sept2017.seuratobject <- CreateSeuratObject(counts = YCC10.DMSO.sept2017.data, project = "YCC10.DMSO.sept2017", min.cells = 3, min.features = 500)

# # Generating quality metrics and selecting (filtering out) cells  -----------------------------
# Add number of genes per UMI for each cell to metadata
YCC10.DMSO.sept2017.seuratobject$log10GenesPerUMI <- log10(YCC10.DMSO.sept2017.seuratobject$nFeature_RNA) / log10(YCC10.DMSO.sept2017.seuratobject$nCount_RNA)

## Calculate percentage of counts originating from the set of features with MT- as set of mitochondrial genes
# exclude cells that are damaged by checking relative expression of mitochondrially derived genes
YCC10.DMSO.sept2017.seuratobject[["percent.mt"]] <- PercentageFeatureSet(YCC10.DMSO.sept2017.seuratobject, pattern = "^MT-")

# Create metadata dataframe as variable to insert addtl metrics for QC without affecting OG merged object
metadata <- YCC10.DMSO.sept2017.seuratobject@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Add sample names
metadata$sample <- NA
metadata$sample<- "YCC10-DMSO-sept2017"

# Add metadata back to Seurat object
YCC10.DMSO.sept2017.seuratobject@meta.data <- metadata

# Add more information to dataset
YCC10.DMSO.sept2017.seuratobject$CellLine <- 'YCC10'
YCC10.DMSO.sept2017.seuratobject$Status <- 'Cancer'
YCC10.DMSO.sept2017.seuratobject$Subtype <- 'EBV'
YCC10.DMSO.sept2017.seuratobject$Treatment <- 'DMSO'
YCC10.DMSO.sept2017.seuratobject$Date <- 'sept2017'
YCC10.DMSO.sept2017.seuratobject$TP53Mutation <- 'No'
YCC10.DMSO.sept2017.seuratobject$Gender <- 'Male'
YCC10.DMSO.sept2017.seuratobject$Culture <- 'Mono'
YCC10.DMSO.sept2017.seuratobject$CellLine.Treatment<- "YCC10_DMSO"
YCC10.DMSO.sept2017.seuratobject$CellLine.Culture <- "YCC10_Mono"
YCC10.DMSO.sept2017.seuratobject$CellLine.Treatment.Culture <- "YCC10_DMSO_Mono"
# ## Assess Quality Metric Plots (save each plot)

# Visualize UMI counts (transcripts) per cell
# Visualize the number UMIs/transcripts per cell
png('nUMIs per cell.png')
metadata %>%
   ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) +
   geom_density(alpha = 0.2) +
   scale_x_log10() +
   theme_classic() +
   ylab("Cell density") +
   geom_vline(xintercept = 6000)+
   ggtitle("nUMIs per cell")
dev.off()

# Visualize the distribution of genes detected per cell via histogram
png('Distn Genes per cell.png')
metadata %>%
 ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) +
 geom_density(alpha = 0.2) +
 theme_classic() +
 scale_x_log10() +
 geom_vline(xintercept = 2000)+
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
png('Before QC nGenes vs nUMIs.png')
metadata %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
  geom_point() +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 8500) +
  geom_hline(yintercept = 2500) +
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
   geom_vline(xintercept = 5.0) +
   ggtitle("Mitochondrial Genes per Cell")
dev.off()

# Visualize the overall novelty of the gene expression by visualizing the genes detected per UMI
png('Novelty plot.png')
metadata %>%
 ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
 geom_density(alpha = 0.2) +
 theme_classic() +
 geom_vline(xintercept = 0.75) +
 ggtitle("Genes detected per UMI")
dev.off()
 # VlnPlot(YCC10.DMSO.sept2017.seuratobject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png('VlnPlot Genes per cell.png')
VlnPlot(YCC10.DMSO.sept2017.seuratobject, features = c("nFeature_RNA") ) + ylab("Expression Level")
dev.off()
png('VlnPlot nCount_RNA per cell.png')
VlnPlot(YCC10.DMSO.sept2017.seuratobject, features = c("nCount_RNA") ) + ylab("Expression Level")
dev.off()
png('VlnPlot percentmt per cell.png')
VlnPlot(YCC10.DMSO.sept2017.seuratobject, features = c("percent.mt") ) + ylab("Expression Level") + scale_y_continuous(breaks=seq(0,100,5))
dev.off()

png('Before QC VlnPlot percentmt per cell.png')
VlnPlot(YCC10.DMSO.sept2017.seuratobject, features = c("percent.mt"),  y.max = 100) +
  ggtitle("Distribution of cells according to the percent of metochondrial reads in cell") +
  ylab("Expression Level") +
  scale_y_continuous(limits=c(0,100)) +
  theme(axis.title = element_text(size=10),axis.text.x = element_text(size=10,angle = 0, hjust = 0.5),axis.text.y = element_text(size=10),plot.title=element_text(size=10,face="bold")) +
  geom_hline(yintercept = 25) +
  theme(legend.position = 'none')
dev.off()

# visualize feature-feature relationships with FeatureScatter
# Save manually
# plot1 <- FeatureScatter(YCC10.DMSO.sept2017.seuratobject, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(YCC10.DMSO.sept2017.seuratobject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))

# Filter out low quality reads using selected thresholds----------------------
YCC10.DMSO.sept2017.seuratobject <- subset(YCC10.DMSO.sept2017.seuratobject,
  subset = nFeature_RNA > 2500 & nFeature_RNA < 8500 &
   percent.mt < 25)
metadata <- YCC10.DMSO.sept2017.seuratobject@meta.data

png('after QC nGenes vs nUMIs.png')
metadata %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
  geom_point() +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  # geom_vline(xintercept = 8500) +
  # geom_hline(yintercept = 2500) +
  facet_wrap(~sample) +
  ggtitle("nGenes vs nUMIs")
dev.off()
png('After QC VlnPlot percentmt per cell.png')
VlnPlot(YCC10.DMSO.sept2017.seuratobject, features = c("percent.mt"),  y.max = 100) +
  ggtitle("Distribution of cells according to the percent of metochondrial reads in cell") +
  ylab("Expression Level") +
  scale_y_continuous(limits=c(0,100)) +
  theme(axis.title = element_text(size=10),axis.text.x = element_text(size=10,angle = 0, hjust = 0.5),axis.text.y = element_text(size=10),plot.title=element_text(size=10,face="bold")) +
  geom_hline(yintercept = 25) +
  theme(legend.position = 'none')
dev.off()
# Run SCTransform: ----------------------
# Normalizes data, Scales Data, Finds variable features
# and regress out percent.mt etc
# SCTransform operate on original unnormalized data in 'RNA' assay slot
#  and saves normalized data in 'SCT' slot
YCC10.DMSO.sept2017.seuratobject <-SCTransform(YCC10.DMSO.sept2017.seuratobject,
  assay = 'RNA',
  new.assay.name = 'SCT',
  vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'), verbose = FALSE)

# Score cells for cell cycle
load("/scrnaseq/data/cycle.rda")
YCC10.DMSO.sept2017.seuratobject <- CellCycleScoring(YCC10.DMSO.sept2017.seuratobject,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes,
                                 assay = 'SCT',
                                 set.ident = TRUE)

# Perform PCA to check cell cycle effect
YCC10.DMSO.sept2017.seuratobject <- RunPCA(YCC10.DMSO.sept2017.seuratobject)

# Plot the PCA colored by cell cycle phase
# Should observe variation by cell cycle
png('Cell Cycle PCA.png')
DimPlot(YCC10.DMSO.sept2017.seuratobject,
       reduction = "pca",
       group.by = "Phase")
dev.off()

# Run SCTransform:
# This time, operate on original unnormalized data in 'RNA' assay slot and
# save to 'SCT' slot, overwriting previous SCT.
YCC10.DMSO.sept2017.seuratobject <-SCTransform(YCC10.DMSO.sept2017.seuratobject,
  assay = 'RNA',
  new.assay.name = 'SCT',
  vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))

# Perform PCA
YCC10.DMSO.sept2017.seuratobject <- RunPCA(YCC10.DMSO.sept2017.seuratobject, features = VariableFeatures(object = YCC10.DMSO.sept2017.seuratobject))

# Now, a PCA on the variable genes no longer returns components associated with cell cycle
png('Cell Cycle Regressed PCA.png')
DimPlot(YCC10.DMSO.sept2017.seuratobject,
        reduction = "pca",
        group.by= "Phase")
dev.off()

# Save SCtransformed object------------------------------------
# To be used as input to integration. RunPCA applied.

save(YCC10.DMSO.sept2017.seuratobject, file="SCT-YCC10-DMSO-sept2017.RData")

# # Linear Dimensional Reduction------------------------------------
pdf('YCC10-DMSO-sept2017-PCAPlots.pdf')
VizDimLoadings(YCC10.DMSO.sept2017.seuratobject, dims = 1:2, reduction = "pca") # all cells plotted on first two PC
DimPlot(YCC10.DMSO.sept2017.seuratobject, reduction = "pca")
DimHeatmap(YCC10.DMSO.sept2017.seuratobject, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(YCC10.DMSO.sept2017.seuratobject, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

# Determine ‘dimensionality’ of the dataset  ----------------------
# The bend at the elbow plot at X PCs, suggests the majority of signal is captured in the first X PCs.
# Choose how many PC dimensions you want to include in the next step of clustering based on the elbow plot.
png('YCC10-DMSO-sept2017-ElbowPlot.png')
ElbowPlot(YCC10.DMSO.sept2017.seuratobject)
dev.off()
