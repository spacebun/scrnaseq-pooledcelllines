This directory contains the scripts that were used to process and analyze the gastric cancer cell lines provided for this project.


First, convert data in `.FASTQ` format to gene expression matrices with scripts in the folder `cellranger`.

The remaining scripts should be called in the following order, in R: 

1. `dataset preparation`
   - For each sample, read in gene expression matrix and perform quality control, SCTRansform, RunPCA.

2. `cell line assignment`
   - For the samples in co-culture, score and predict cells for similarity to SNU719 and NCC24 using the data of SNU719 and NCC24 in monoculture. Assign SNU719 and NCC24 cells in co-culture as well.
3. `analysis`
   - `verify cell line assigment`: Use similarity in pattern of differential gene expression between SNU719 and NCC24 in monoculture to SNU719 and NCC24 in co-culture, to assess the accuracy of cell line assignment. 
   - Perform 2-way differential expression analysis of SNU719 cells across treatment condition or culture condition. Obtain list of differentially expressed genes and the heatmap of the top 40 most differentially expressed genes.


