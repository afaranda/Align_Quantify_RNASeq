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
# TO DO:
#      -  Determine whether to use the strand "+/-" value from the first or
#         or the second read in a pair
#      -  Create and run a mock vaildation set
#
###############################################################################


## Setup Environment
export ALIGNDIR=$(pwd)/Alignments
export COUNTDIR=$(pwd)/Counts
export BEDPATH=/work/abf/MouseEnsembl101/rseqc_gene_models.bed
export GTFPATH=/work/abf/MouseEnsembl101/Mus_musculus.GRCm38.101.gtf
export fasta=/work/abf/MouseEnsembl101/Mus_musculus.GRCm38.dna.primary_assembly.fa
export MINQUAL=30
export MINGC=0.7

## Iterate Over Alignments; generate and index name sorted BAM files
for b in $(find $ALIGNDIR -type f -name "*sorted_alignment.bam")
do
    fn=${b/_sorted_alignment\.bam/_byname_alignment.bam}
    if [ ! -f $fn ]
    then
       samtools sort -n -@ 4 -m 12G ${b} -o $fn
    fi
done

## Extract Ranges from name sorted BAM files in bedpe format
## convert bedpe to BED spanning full fragment, calculate GC for each fragment
## and identify genes that overlap each fragment

for b in $(find $ALIGNDIR -type f -name "*byname_alignment.bam")
do
    obed=$(echo $b | sed 's/byname_alignment\.bam/high_gc_frags.bed/')
    ofn=$(echo $b | sed 's/byname_alignment\.bam/high_gc_frags_by_gene.txt/')
    oht=$(echo $b | sed 's/byname_alignment\.bam/high_gc_frags_htseq_count.txt/')
    ID=${b/_byname_alignment.bam/}
    ID=${ID##*/}
    if [ ! -f $obed ]
    then
	echo input bam file: $b
	echo output bed file: $obed
	bedtools bamtobed -bedpe -i <(samtools view -q $MINQUAL -h -f3 -F 256 $b)\
	    | gawk -v FS="\t"\
		   -v OFS="\t"\
                   '\
                   ( $1 !~ /\./ && $6 > $2)\
                   { 
	               print $1, $2, $6, $7, $8, $10
                   }
                   '\
	    | bedtools nuc -s -fi $fasta -bed stdin\
	    | gawk -v OFS="\t" -v MG=$MINGC '($8 < MG) {print $1, $2, $3, $4, $8, $6}'\
	    | tee -a ${obed}\
	    | bedtools intersect -s -c -b stdin -a <(\
	       gawk -v OFS="\t"\
	       	    -v FS="\t|; "\
	       	    '{gsub("gene_id ","",$9); gsub("\042","",$9)}
                     {gsub("gene_name ","",$11); gsub("\042","",$11)}
                     ($3 ==  "gene")\
                     {print $1, $2, $3, $4, $5, $6, $7, $8,$9"\t", $11}'\
	       	     ${GTFPATH}) > ${ofn}
    else
	echo Already extracted high gc bed file for ${ID}
    fi
    
    obam=${obed/high_gc_frags.bed/sorted_highGC_alignment.bam}
    if [ ! -f $obam ]
    then
       ibam=${b/byname_alignment\.bam/sorted_alignment.bam}
       echo Filtering position-sorted alignment for $ID
       echo $ibam
       echo $obed
       echo $obam
       samtools view -h -b -L $obed $ibam > $obam
       samtools index $obam
    else
	echo Already generated high-gc filtered alignments for ${ID}
    fi

    
    if [ ! -f ${COUNTDIR}/${ID}_HighGC_Count.txt ]
    then
	htseq-count \
	    -i gene_id -r pos -f bam -s reverse -m union --type exon \
	    $obam \
	    $GTFPATH > ${COUNTDIR}/${ID}_HighGC_Count.txt
    else
	echo Already counted high gc reads with htseq_count for ${ID}
    fi
    echo
done
