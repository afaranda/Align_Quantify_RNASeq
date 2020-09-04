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
export ALIGNDIR=../30J*/Alignments
export BEDPATH=/work/abf/MouseEnsembl100/rseqc_gene_models.bed
export GTFPATH=/work/abf/MouseEnsembl100/Mus_musculus.GRCm38.100.gtf
export fasta=/work/abf/MouseEnsembl100/Mus_musculus.GRCm38.dna.primary_assembly.fa
export MINQUAL=30
export MINGC=0.7

## Iterate Over Alignments; generate and index name sorted BAM files
for b in $(find $ALIGNDIR -type f -name "*sorted_alignment.bam")
do
    fn=${b/_sorted_alignment\.bam/_byname_alignment.bam}
    if [ ! -f $fn ]
       samtools sort -n -@ 4 -m 12G ${b} -o $fn
    fi
done

## Extract Ranges from name sorted BAM files in bedpe format
## convert bedpe to BED spanning full fragment, calculate GC for each fragment
## and identify genes that overlap each fragment
for b in$(find $ALIGNDIR -type f -name "*byname_alignment.bam")
do
    bedtools bamtobed -bedpe -i $b\
	| gawk -v FS="\t"\
	       -v OFS="\t"\
	       -v MQ=$MINQUAL\
        '\
        ( $1 !~ /\./ && $9 >= MQ)\
        {
             print $1, $2, $6, $7, $8, $9
        }
        '\
	| bedtools nuc -fi $fasta -i stdin
done
