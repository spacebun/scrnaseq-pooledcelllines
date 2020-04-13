#!/bin/bash/
# This script takes a fastq file of scRNA-Seq data,
# and aligns sequencing reads in FASTQ files to a reference transcriptome
# and generates a .cloupe file for visualization and analysis in Loupe Cell Browser,
# along with a number of other outputs compatible with other publicly-available tools for further analysis.

# USAGE: qsub -b n -pe smp 8 /gpfs/serena/scrnaseq/scripts/cellranger/mixed/cellranger_count_mixed_Aza.sh

# output folder: /gpfs/serena/scrnaseq/results/run_cellranger_count_mixed/mixed_Aza
# data directory: /gpfs/serena/mixed_cell_lines/cell_Aza
# sample id: cell_Aza

#go to directory to run the analysis in
cd /gpfs/serena/scrnaseq/results/run_cellranger_count_mixed

#run cellranger count
# --id=<name of output folder created by cell ranger count. Recommended to name as desired sample name.>
# --fastqs=<data directory>
# --sample=<specify which samples to use based off of the sample id at the beginning of the FASTQ file name
# --transcriptome=path to the transcriptome reference package
cellranger count --id=mixed_Aza \
--fastqs=/gpfs/serena/mixed_cell_lines/cell_Aza \
--sample=cell_Aza \
--transcriptome=/gpfs/serena/apps/cellranger_references/refdata-cellranger-GRCh38-3.0.0
