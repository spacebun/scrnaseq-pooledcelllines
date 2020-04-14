This directory contains the scripts for the creation of pre-processed Seurat
objects for the cell lines in co-culture, treated with DMSO and Aza as listed below.

The subdirectories are organised as follows:

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
