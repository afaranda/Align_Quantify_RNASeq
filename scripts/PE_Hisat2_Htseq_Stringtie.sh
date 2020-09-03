#!/bin/bash
#SBATCH --job-name=PE_Hisat2_Htseq_Stringtie
#SBATCH --mem=128000
#SBATCH --ntasks=8

ID=$(echo $1 | sed 's/_L[0-9]\{3\}_R[12]_[0-9]\{3\}\.fastq\.gz//')
F1=$(echo $1 | sed 's/\.fastq\.gz/_fastqc.html/')
echo $ID

## Pretrim FastQC
if [ -f ${PRETRIM_QC}/${F1} ]
then
    echo skipping pretrim FastQC for $1
    echo skipping pretrim FastQC for $2
else
    fastqc ${FASTQDIR}/${1} -o ${PRETRIM_QC}
    fastqc ${FASTQDIR}/${2} -o ${PRETRIM_QC}
fi

## Get Trimmed Filenames
R1=$(echo $1 | sed 's/\.fastq\.gz/_val_1.fq.gz/')
R2=$(echo $2 | sed 's/\.fastq\.gz/_val_2.fq.gz/')

## Trimgalore
if [ -f ${TRIMDIR}/${R1} ]
then
    echo skipping trim for $ID
else
    trim_galore -o $TRIMDIR \
		--length 85 \
		--paired \
		--clip_R2 3 \
		--phred33 \
		${FASTQDIR}/${1} \
		${FASTQDIR}/${2}
fi

## Run Posttrim FastQC
R1=$(echo $1 | sed 's/\.fastq\.gz/_val_1.fq.gz/')
R2=$(echo $2 | sed 's/\.fastq\.gz/_val_2.fq.gz/')


fastqc ${TRIMDIR}/${R1} -o ${POSTTRIM_QC}
fastqc ${TRIMDIR}/${R2} -o ${POSTTRIM_QC}

## Align Trimmed Reads Using Hisat2
hisat2 -p8\
       --verbose\
       --phred33\
       --dta\
       --fr\
       -x $HISAT2_PREFIX\
       -1 ${TRIMDIR}/${R1}\
       -2 ${TRIMDIR}/${R2}\
       -S ${ALIGNDIR}/${ID}_aligned_reads.sam


## compress, sort, and index alignments
samtools view -bS ${ALIGNDIR}/${ID}_aligned_reads.sam > ${ALIGNDIR}/${ID}_aligned_reads.bam

samtools sort -o ${ALIGNDIR}/${ID}_sorted_alignment.bam ${ALIGNDIR}/${ID}_aligned_reads.bam

samtools index ${ALIGNDIR}/${ID}_sorted_alignment.bam

rm ${ALIGNDIR}/${ID}_aligned_reads.bam ${ALIGNDIR}/${ID}_aligned_reads.sam

## Get Gene-level counts with Htseq Count
# Call htseq-count on the target bam file
htseq-count \
     -i gene_id -r pos -f bam -s reverse -m union --type exon \
     ${ALIGNDIR}/${ID}_sorted_alignment.bam \
     $GTFPATH > ${COUNTDIR}/${ID}_GeneCount.txt


# Call stringtie on the target bam file
stringtie ${ALIGNDIR}/${ID}_sorted_alignment.bam -p8\
     --fr\
     -e -B\
     -G $GTFPATH\
     -o ${STRNGDIR}/${ID}/${ID}.gtf\
     -A ${STRNGDIR}/${ID}/${ID}.txt
