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
export INDEXDIR=                    # Path to Genome Index Directory
export ALIGNDIR=$(pwd)/Alignments   # Path to alignment output directory
export COUNTDIR=$(pwd)/Counts       # Path to count output directory
export STRNGDIR=$(pwd)/Stringtie    # Path to Stringtie output directory
export QCDIR=$(pwd)/QC_Results      # Path to QC Output Directory
export POSTTRIM_QC=$(pwd)/PT_FastQC # Path to Post Trimming FastQC(unify later)
export GTFPATH=    # Path to GTF File
export TRIMDIR=$(pwd)/Trimmed   # Path to Trimmed Reads Directory
export FQTARGET="L[0-9]\{3\}_R1_[0-9]\{3\}\.fastq\.gz"  # Regex for fastq files

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
	 $QCDIR \
	 $POSTTRIM_QC \
	 $TRIMDIR)
for D in ${ARR[@]}; do
    if [ ! -d $D ]; then
	mkdir $D
	
    else
	rm -rf $D
	mkdir $D
    fi
done
    
# Iterate over FastQ Files, launch aprocessing script for each pair of reads
echo $(ls $FASTQDIR)
for R1 in $(ls $FASTQDIR | grep $FQTARGET)
do
    CK=$(echo $R1 | sed 's/\.fastq\.gz//')
    echo "Running" ${CK}
    if [ ! -f ${POSTTRIM_QC}/${CK} ] ; then
	R2=$(echo $R1 | sed 's/R1/R2/')
	echo "Running" sbatch PE_Hisat2_Htseq_Stringtie.sh $R1 $R2
	./scripts/PE_Hisat2_Htseq_Stringtie.sh $R1 $R2
    fi
done

# Aggregate QC Results in One Place
