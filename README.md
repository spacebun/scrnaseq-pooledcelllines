# Final Year Project: Bioinformatics for Virtual Sorting of Cancer Cells in Mixed Populations

## Description
This repository contains the scripts and results of sorting and analysing scRNA-seq of gastric cancer cell lines in monoculture and in a pooled co-culture. 

### Scripts
The sub-folders in https://github.com/spacebun/scrnaseq-pooledcelllines/tree/master/scripts represents the phases of the pipeline. 

| s/n | Phase | Purpose |
| -----| ----- | -------| 
| 1 | dataset preparation | For each sample, read in gene expression matrix and perform quality control, SCTRansform, RunPCA.|
| 2 | cell line assignment | For the samples in co-culture, score and predict cells for similarity to SNU719 and NCC24 using the data of SNU719 and NCC24 in monoculture. Assign SNU719 and NCC24 cells in co-culture as well. |
| 3 | analysis | 1. verify cell line assigment: Use similarity in pattern of differential gene expression between SNU719 and NCC24 in monoculture to SNU719 and NCC24 in co-culture, to assess the accuracy of cell line assignment.      2. Perform 2-way differential expression analysis of SNU719 cells across treatment condition or culture condition. Obtain list of differentially expressed genes and the heatmap of the top 40 most differentially expressed genes. |



#### Packages for R

**(1)** The following packages listed below are from **CRAN** and were installed using the `install.packages()` function. 

1. `Seurat`
1. `sctransform`
1. `tidyverse`
1. `dplyr`
1. `ggplot2`
1. `ggrepel`
1. `cowplot`
1. `Matrix`
1. `scales`
1. `bitops`
1. `RCurl`
1. `Matrix.utils`
1. `RColorBrewer`


**(2)** The following package from **Bioconductor** was installed using the the `BiocManager::install()` function.

1. `SingleCellExperiment`

