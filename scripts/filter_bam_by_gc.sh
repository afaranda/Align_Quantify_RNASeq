#!/bin/bash
#SBATCH --job-name=filter_bam_by_gc
#SBATCH --mem=64000
#SBATCH --ntasks=4

###############################################################################
# 
# File: filter_bam_by_gc.sh
# Purpose: Select reads from a BAM file based on their GC content and
#          and identify the genes / features that they overlap. 
# Created: August 26, 2020
# Author: Adam Faranda
#
###############################################################################


## Setup Environment
export ALIGNDIR=$(pwd)/Alignments
export BEDPATH=/work/abf/MouseEnsembl100/rseqc_gene_models.bed
export GTFPATH=/work/abf/MouseEnsembl100/Mus_musculus.GRCm38.100.gtf

## Iterate Over Alignments; generate and index name sorted BAM files
for b in $(find ${ALIGNDIR} -type f -name "*sorted_alignment.bam")
do
samtools sort 
done


## Extract Ranges from name sorted BAM files in bedpe format
## convert bedpe to BED spanning full fragment, calculate GC for each fragment
## and identify genes that overlap each fragment
