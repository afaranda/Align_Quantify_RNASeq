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
F1=$(echo $R1 | sed 's/\.fq\.gz/_fastqc.html/')
if [ -f ${POSTTRIM_QC}/${F1} ]
then
    echo skipping posttrim QC $R1
    echo skipping posttrim QC $R2
else
    fastqc ${TRIMDIR}/${R1} -o ${POSTTRIM_QC}
    fastqc ${TRIMDIR}/${R2} -o ${POSTTRIM_QC}
fi

## Align Trimmed Reads Using Hisat2
if [ -f ${ALIGNDIR}/${ID}_sorted_alignment.bam ]
then
   echo skipping alignment for sample $ID
else
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
fi

## Get Gene-level counts with Htseq Count
# Call htseq-count on the target bam file
if [ -f ${COUNTDIR}/${ID}_GeneCount.txt ]
then
    echo skipping htseq-count for $ID
else
    htseq-count \
	-i gene_id -r pos -f bam -s reverse -m union --type exon \
	${ALIGNDIR}/${ID}_sorted_alignment.bam \
	$GTFPATH > ${COUNTDIR}/${ID}_GeneCount.txt
fi

# Call stringtie on the target bam file
if [ -f ${STRTIEDIR}/${ID}/${ID}.txt ]
then
    echo skipping stringtie for $ID
else
    stringtie ${ALIGNDIR}/${ID}_sorted_alignment.bam -p8\
	      --fr\
	      -e -B\
	      -G $GTFPATH\
	      -o ${STRNGDIR}/${ID}/${ID}.gtf\
	      -A ${STRNGDIR}/${ID}/${ID}.txt
fi
