This directory contains the scripts for the creation of pre-processed Seurat
objects for each of the pure cell lines as listed below.

The subdirectories are organised as follows:

* april_2017
⋅⋅⋅SNU719-Aza-april2017
⋅⋅⋅You can have properly indented paragraphs within list items. Notice the blank line above, and the leading spaces (at least one, but we'll use three here to also align the raw Markdown).

⋅⋅⋅To have a line break without a paragraph, you will need to use two trailing spaces.⋅⋅
⋅⋅⋅Note that this line is separate, but within the same paragraph.⋅⋅
⋅⋅⋅(This is contrary to the typical GFM line break behaviour, where trailing spaces are not required.)



april_2017
	SNU719-Aza-april2017
  SNU719-DMSO-april2017

july_2017
	HFE145-july2017
	SNU719-DMSO-july2017

sept_2017
	NCC24-DMSO-sept2017
	SNU719-Aza-sept2017
	YCC10-Aza-sept2017
	YCC10-DMSO-sept2017


Overview of each script:



Steps:
Read in data and create Seurat object
Create, plot and assess quality metrics of data (mainly percent.mt)
Perform SCTransform, vars.to.regess = c("percent.mt, Filter out low quality cells/reads from Seurat object
Perform CellCycleScoring, FindVariableFeatures, ScaleData, RunPCA, then plot to determine if cell cycle is major source of variation
(include these 2 lines of code: marrow <- RunPCA(object, features = c(s.genes, g2m.genes))
DimPlot(object)
)
Always regress out variation due to cell cycle with
Perform SCTransform, vars.to.regess = c("percent.mt, Filter out low quality cells/reads from Seurat object
