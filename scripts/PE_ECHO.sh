#!/bin/bash
#SBATCH --job-name=PE_Hisat2_Htseq_Stringtie
#SBATCH --mem=128000
#SBATCH --ntasks=8

ID=$(echo $1 | sed 's/_L[0-9]\{3\}_R[12]_[0-9]\{3\}\.fastq\.gz//')
echo SAMPLE $ID

## Pretrim FastQC
echo fastqc ${FASTQDIR}/${1} -o ${QCDIR}
echo fastqc ${FASTQDIR}/${2} -o ${QCDIR}

## Trimgalore
echo trim_galore -o $TRIMMED \
	    --length 85 \       # Ajust based on number of reads
	    --paired \
	    --clip_R2 3 \
	    --phred33 \
	    ${FASTQDIR}/${1} \
	    ${FASTQDIR}/${2}
	    

## Run Posttrim FastQC
R1=$(echo $1 | sed 's/\.fastq\.gz/_val_1.fq.gz/')
R2=$(echo $2 | sed 's/\.fastq\.gz/_val_2.fq.gz/')


echo fastqc ${TRIMDIR}/${R1} -o ${POSTTRIM_QC}
echo fastqc ${TRIMDIR}/${R2} -o ${POSTTRIM_QC}

## Align Trimmed Reads Using Hisat2
echo "hisat2 -p8\
       --verbose\
       --phred33\
       --dta\
       --rf\                   # Adjust based on library method
       -x genome_tran\                   
       -1 ${FASTQDIR}/${R1}\
       -2 ${FASTQDIR}/${R2}\
       -S ${ALIGNDIR}/${ID}_aligned_reads.sam"


## compress, sort, and index alignments
echo "samtools view -bS ${ALIGNDIR}/${ID}_aligned_reads.sam > ${ALIGNDIR}/${ID}_aligned_reads.bam"

echo samtools sort -o ${ALIGNDIR}/${ID}_sorted_alignment.bam ${ALIGNDIR}/${ID}_aligned_reads.bam

echo samtools index ${ALIGNDIR}/${ID}_sorted_alignment.bam

rm ${ALIGNDIR}/${ID}_aligned_reads.bam ${ALIGNDIR}/${ID}_aligned_reads.sam

## Get Gene-level counts with Htseq Count
# Call htseq-count on the target bam file
echo "htseq-count -i gene_id -f bam -s reverse -m union --type exon $1 $GTF_FILE > $2"
echo "htseq-count \
     -i gene_id -r pos -f bam -s reverse -m union --type exon \
     ${ALIGNDIR}/${ID}_sorted_alignment.bam \
     $GTFPATH > ${COUNTDIR}/${ID}_GeneCount.txt"


# Call stringtie on the target bam file

echo "stringtie $ID -p8\
     --rf\                               # Adjust based on library method
     -e -B\
     -G $GTFPATH\
     -o ${STRTIEDIR}/${ID}/${ID}.gtf\
     -A ${STRTIEDIR}/${ID}/${ID}.txt"

