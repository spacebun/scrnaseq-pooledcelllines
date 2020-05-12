This directory contains the results of analysis of the gastric cancer cell lines provided in monoculture and co-culture.

The subdirectories are organised as follows:

* `dataset preparation`: For each sample, read in gene expression matrix and perform quality control, SCTRansform, RunPCA.
  * `mono`
    * `april_2017`
      * SNU719-Aza-april2017
      * SNU719-DMSO-april2017
    * `july_2017`
      * HFE145-july2017
      * SNU719-DMSO-july2017
    * `sept_2017`
      * NCC24-DMSO-sept2017
      * SNU719-Aza-sept2017
      * YCC10-Aza-sept2017
      * YCC10-DMSO-sept2017
  * co
    * co-DMSO
    * co-Aza
    
* `cell line assignment`: For the samples in co-culture, score and predict cells for similarity to SNU719 and NCC24 using the data of SNU719 and NCC24 in monoculture. Assign SNU719 and NCC24 cells in co-culture as well.
  * `DMSO`: Assignment of cells in co-culture treated with DMSO.
  * `Aza`Assignment of cells in co-culture treated with Aza.

*  `analysis`: Contains heatmaps and list of differentially expressed genes from doing two-group differential expression testing of the samples.
  * `verify cell line assignment`: Comparing SNU719 and NCC24 in monoculture. Then, use the top 20 differentially expressed genes from this comparison to assess similarity in pattern of differential gene expression pattern of SNU719 and NCC24 in co-culture to that of SNU719 and NCC24 in monoculture.
  * `monoDMSOvsmonoAza`: Comparing the Aza-treated and DMSO-treated SNU719 cells in monoculture.
  * `coDMSOvscoAza`: Comparing the Aza-treated and DMSO-treated SNU719 cells in co-culture.
  * `monoDMSOvscoDMSO`: Comparing the DMSO-treated SNU719 cells in monoculture and in co-culture.
  * `monoAzavscoAza`: Comparing the Aza-treated SNU719 cells in monoculture and in co-culture.


