#!/bin/bash/

#This is a version done on 16 Dec as an error was found in the outs of the run of the original script.

# This script takes a fastq file of RNA-Seq data,
# and aligns sequencing reads in FASTQ files to a reference transcriptome 
# and generates a .cloupe file for visualization and analysis in Loupe Cell Browser, 
# along with a number of other outputs compatible with other publicly-available tools for further analysis. 

# USAGE: qsub -b n -pe smp 8 /gpfs/serena/scrnaseq/scripts/cellranger/cellranger_count_SNU719-DMSO-april2017-v2.sh


# data: /gpfs/serena/single_cell/april_2017/sample1-ut
# sample name: sample1-ut-cancer

#go to directory to run the analysis in 
cd /gpfs/serena/scrnaseq/results/run_cellranger_count


cellranger count --id=run_count_SNU719-DMSO-april2017 \
--fastqs=/gpfs/serena/single_cell/april_2017/sample1-ut \
--sample=sample1-ut-cancer \
--transcriptome=/gpfs/serena/apps/cellranger_references/refdata-cellranger-GRCh38-3.0.0