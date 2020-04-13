#!/bin/bash/
# This script takes a fastq file of RNA-Seq data,
# and aligns sequencing reads in FASTQ files to a reference transcriptome 
# and generates a .cloupe file for visualization and analysis in Loupe Cell Browser, 
# along with a number of other outputs compatible with other publicly-available tools for further analysis. 

# USAGE: qsub -V -b n -cwd /gpfs/serena/rnaseq_test/cellranger_count_NCC24-DMSO-sept2017.sh

# output folder: run_count_NCC24-DMSO-sept2017
# data directory: /gpfs/serena/single_cell/sept_2017/10X_K6_NCC24
# sample id: 10X_K6_NCC24
 
#go to directory to run the analysis in
#mkdir /gpfs/serena/rnaseq_test/run_cellranger_count
cd /gpfs/serena/rnaseq_test/run_cellranger_count

#run cellranger count 
# --id=<name of output folder created by cell ranger count>
# --fastqs=<data directory>
# --sample=<specify which samples to use based off of the sample id at the beginning of the FASTQ file name
# --transcriptome=path to the transcriptome reference package
cellranger count --id=run_count_NCC24-DMSO-sept2017 \
--fastqs=/gpfs/serena/single_cell/sept_2017/10X_K6_NCC24 \
--sample=10X_K6_NCC24 \
--transcriptome=/gpfs/serena/apps/cellranger_references/refdata-cellranger-GRCh38-3.0.0