Contains scripts that read in gene expression matrices and performs quality control, normalization, scaling, and PCA. 

This directory contains the scripts for the creation of pre-processed Seurat
objects for each of the pure cell lines as listed below.

The subdirectories are organised as follows:
* pure_indiv
  * april_2017
    * SNU719-Aza-april2017
    * SNU719-DMSO-april2017

  * july_2017
    * HFE145-july2017
    * SNU719-DMSO-july2017

  * july_2017
    * NCC24-DMSO-sept2017
    * SNU719-Aza-sept2017
    * YCC10-Aza-sept2017
    * YCC10-DMSO-sept2017
    
* mixed_indiv
  * mixed_DMSO
  * mixed_Aza
  
Overview of each script:
1. Read in data and create Seurat object, adding information about sample
2. Create, plot and assess quality metrics of data 
3. Filter out poor quality cells
4. Perform SCTransform
5. Perform CellCycleScoring, FindVariableFeatures, ScaleData, RunPCA, then plot to determine if cell cycle is major source of variation
6. Perform SCTransform, regressing out variation due to cell cycle 
7. RunPCA on object
8. Save object
