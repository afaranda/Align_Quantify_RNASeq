#!/bin/bash
#SBATCH --job-name=PE_Hisat2_Htseq_Stringtie
#SBATCH --mem=128000
#SBATCH --ntasks=8


## Get Trimmed Filenames and sample label
ID=$(echo $1 | sed 's/_L[0-9]\{3\}_R[12]_[0-9]\{3\}\.fastq\.gz//')
R1=$(echo $1 | sed 's/\.fastq\.gz/_val_1.fq.gz/')
R2=$(echo $2 | sed 's/\.fastq\.gz/_val_2.fq.gz/')

## Align Trimmed Reads To Ribosomal Index Using Hisat2
if [ -f ${ALIGNDIR}/${ID}_sorted_ribo.bam ]
then
   echo skipping alignment for sample $ID
else   
    hisat2 -p8\
	   --verbose\
	   --phred33\
	   --dta\
	   --fr\
	   -x Rn45s_Index\
	   -1 ${TRIMDIR}/${R1}\
	   -2 ${TRIMDIR}/${R2}\
	   -S ${ALIGNDIR}/${ID}_ribo_reads.sam

    ## compress, sort, and index alignments
    samtools view -bS ${ALIGNDIR}/${ID}_ribo_reads.sam > ${ALIGNDIR}/${ID}_ribo_reads.bam
    
    samtools sort -o ${ALIGNDIR}/${ID}_sorted_ribo.bam ${ALIGNDIR}/${ID}_ribo_reads.bam
    
    samtools index ${ALIGNDIR}/${ID}_sorted_ribo.bam
    
    rm ${ALIGNDIR}/${ID}_ribo_reads.bam ${ALIGNDIR}/${ID}_ribo_reads.sam
fi

## 

