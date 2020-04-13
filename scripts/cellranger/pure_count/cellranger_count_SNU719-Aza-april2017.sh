#!/bin/bash/
# This script takes a fastq file of RNA-Seq data,
# and aligns sequencing reads in FASTQ files to a reference transcriptome 
# and generates a .cloupe file for visualization and analysis in Loupe Cell Browser, 
# along with a number of other outputs compatible with other publicly-available tools for further analysis. 

# USAGE: qsub -V -b n -cwd /gpfs/serena/rnaseq_test/cellranger_count_SNU719-Aza-april2017.sh

# data: /gpfs/serena/single_cell/april_2017/sample1-t
# sample name: sample1-t-cancer
 
#making a directory to run the analysis in
mkdir /gpfs/serena/rnaseq_test/run_cellranger_count
cd /gpfs/serena/rnaseq_test/run_cellranger_count


cellranger count --id=run_count_SNU719-Aza-april2017 \
--fastqs=/gpfs/serena/single_cell/april_2017/sample1-t \
--sample=sample1-t-cancer \
--transcriptome=/gpfs/serena/apps/cellranger_references/refdata-cellranger-GRCh38-3.0.0