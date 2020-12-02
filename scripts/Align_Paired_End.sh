#!/bin/bash
#
# File: Process_Paired_End.sh
# Author: Adam Faranda
# Purpose:
#     Align paired end RNASeq reads to a refference genome and quantify
#     gene / transcript abundance
# Created: July 27, 2020



# Define Main Directories
export PATH=${PATH}:$(pwd)/scripts  # Add this pipeline to the executeable path
export FASTQDIR=$(pwd)/fastq        # Path to directory with reads
export HISAT2_INDEXES=/work/abf/MouseEnsembl101                    # Path to Genome Index Directory
export HISAT2_PREFIX=EnsMm_101      # Prefix for Hisat2 index
export ALIGNDIR=$(pwd)/Alignments   # Path to alignment output directory
export COUNTDIR=$(pwd)/Counts       # Path to count output directory
export STRNGDIR=$(pwd)/Stringtie    # Path to Stringtie output directory
export PRETRIM_QC=$(pwd)/PRE_FastQC      # Path to Pre Trim QC Output Directory
export POSTTRIM_QC=$(pwd)/POST_FastQC # Path to Post Trimming FastQC(unify later)
export GTFPATH=/work/abf/MouseEnsembl101/Mus_musculus.GRCm38.101.gtf    # Path to GTF File
export RSEQCDIR=$(pwd)/RSeQC_Results        # Path to RSeQC results
export BEDPATH=/work/abf/MouseEnsembl101/rseqc_gene_models.bed # Path to bed file for RSeQC
export TRIMDIR=$(pwd)/Trimmed   # Path to Trimmed Reads Directory
export FQTARGET="L[0-9]\{3\}_R1_[0-9]\{3\}\.fastq\.gz"  # Regex for fastq files
export DELOLD=0 # Set to 1 to delete previous results

# De-Gitify and Create directories if none exist
if [ -d .git ]; then
    rm -rf .git
fi

if [ -f .gitignore ]; then
   rm .gitignore
fi

ARR=($ALIGNDIR \
	 $COUNTDIR \
	 $STRNGDIR \
	 $PRETRIM_QC \
	 $POSTTRIM_QC \
	 $TRIMDIR \
	 $RSEQCDIR)
for D in ${ARR[@]}; do
    if [ ! -d $D ]; then
	mkdir $D
	
    elif [ $DELOLD == 1 ]; then
	rm -rf $D
	mkdir $D
    fi
done
    
# Iterate over FastQ Files, launch aprocessing script for each pair of reads
echo $(ls -lh $FASTQDIR)
JOBS=""
for R1 in $(ls $FASTQDIR | grep $FQTARGET)
do
    # Align Paired and reads to genome and quantify abundance
    R2=$(echo $R1 | sed 's/R1/R2/')
    echo "Running" sbatch PE_Hisat2_Htseq_Stringtie.sh $R1 $R2
    JB=$(sbatch PE_Hisat2_Htseq_Stringtie.sh $R1 $R2 | gawk '{print $4}')
    JOBS=${JOBS},afterok:${JB}

    # Align Paired end reads to Ribosomal DNA fragment Rn45s
    RB=$(sbatch\
	     --dependency=afterok:${JB}\
	     PE_ribosomal_content_analysis.sh $R1 $R2 | gawk '{print $4}'
      )
    ROBS=${ROBS},afterok:${RB}
done

# Aggregate QC Results in One Place
export JOBS=$(echo $JOBS | sed 's/,//')
export ROBS=$(echo $ROBS | sed 's/,//')
sbatch --dependency=${JOBS},${ROBS} RunQC.sh
