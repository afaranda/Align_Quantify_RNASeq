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
		--cores 8 \
		${FASTQDIR}/${1} \
		${FASTQDIR}/${2}
fi


## Align Trimmed Reads To Ribosomal Index Using Hisat2
if [ -f ${ALIGNDIR}/${ID}_byname_ribo.bam ]
then
   echo skipping alignment for sample $ID
else   
    hisat2 -p8\
	   --verbose\
	   --phred33\
	   --dta\
	   --fr\
	   --summary-file ${ALIGNDIR}/${ID}_RiboAlignStat.txt\
	   -x Rn45s_Index\
	   -1 ${TRIMDIR}/${R1}\
	   -2 ${TRIMDIR}/${R2}\
	   -S ${ALIGNDIR}/${ID}_ribo_reads.sam

    ## compress, sort, and index alignments
    samtools view -@ 8 -bS ${ALIGNDIR}/${ID}_ribo_reads.sam > ${ALIGNDIR}/${ID}_ribo_reads.bam
    
    samtools sort -n -@ 8 -m 12G -o ${ALIGNDIR}/${ID}_byname_ribo.bam ${ALIGNDIR}/${ID}_ribo_reads.bam
    
    rm ${ALIGNDIR}/${ID}_ribo_reads.bam ${ALIGNDIR}/${ID}_ribo_reads.sam
fi

## Extract reads that Failed to align to the ribosome into
## new fastq files

if [ -f ${TRIMDIR}/${ID}_ribotrim_R1.fq ]
then
    echo skipping alignment for sample $ID
else
    bedtools bamtofastq\
	     -i <(samtools view -h -f 13 ${ALIGNDIR}/${ID}_byname_ribo.bam)\
	     -fq ${TRIMDIR}/${ID}_ribotrim_R1.fq\
	     -fq2 ${TRIMDIR}/${ID}_ribotrim_R2.fq
    
fi

## Run Posttrim FastQC
F1=$(echo $R1 | sed 's/_val_1\.fq\.gz/_ribotrim_R1_fastqc.html/')
if [ -f ${POSTTRIM_QC}/${F1} ]
then
    echo skipping posttrim QC ${ID}_ribotrim_R1.fq
    echo skipping posttrim QC ${ID}_ribotrim_R2.fq
else
    fastqc ${TRIMDIR}/${ID}_ribotrim_R1.fq -o ${POSTTRIM_QC}
    fastqc ${TRIMDIR}/${ID}_ribotrim_R2.fq -o ${POSTTRIM_QC}
fi

## Align Trimmed Reads Using Hisat2
if [ -f ${ALIGNDIR}/${ID}_sorted_alignment.bam ]
then
    echo skipping alignment for sample $ID
else
    echo aligning to genome
    echo read1 file ${ID}_ribotrim_R1.fq
    echo read2 file ${ID}_ribotrim_R2.fq
    
    hisat2 -p8\
	   --verbose\
	   --phred33\
	   --dta\
	   --fr\
	   --summary-file ${ALIGNDIR}/${ID}_AlignStat.txt\
	   -x $HISAT2_PREFIX\
	   -1 ${TRIMDIR}/${ID}_ribotrim_R1.fq\
	   -2 ${TRIMDIR}/${ID}_ribotrim_R2.fq\
	   --rna-strandness RF\
	   -S ${ALIGNDIR}/${ID}_aligned_reads.sam
    
    ## compress, sort, and index alignments
    samtools view -@ 8 -bS ${ALIGNDIR}/${ID}_aligned_reads.sam > ${ALIGNDIR}/${ID}_aligned_reads.bam
    
    samtools sort -@ 8 -m 12G -o ${ALIGNDIR}/${ID}_sorted_rf_alignment.bam ${ALIGNDIR}/${ID}_aligned_reads.bam
    
    samtools index ${ALIGNDIR}/${ID}_sorted_rf_alignment.bam
    
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
	${ALIGNDIR}/${ID}_sorted_rf_alignment.bam \
	$GTFPATH > ${COUNTDIR}/${ID}_rf_GeneCount.txt
fi

## Estimate Transcript Abundances with Kallisto
# Call stringtie on the target bam file
if [ -f ${KLSTODIR}/${ID}/abundances.txt ]
then
    echo skipping stringtie for $ID
else
    mkdir ${KLSTODIR}/${ID}
    kallisto quant -i ${KLSTOIDX} -t 8 -b 25 --seed 8253\
	     -o ${KLSTODIR}/${ID}\
	     ${TRIMDIR}/${ID}_ribotrim_R1.fq\
	     ${TRIMDIR}/${ID}_ribotrim_R2.fq
fi

# Estimate Transcript Integrity using iteratively generated script
tin.py -i ${ALIGNDIR}/${ID}_sorted_rf_alignment.bam -r $BEDPATH
