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
else
    fastqc ${FASTQDIR}/${1} -o ${PRETRIM_QC}
fi

## Get Trimmed Filename
R1=$(echo $1 | sed 's/\.fastq\.gz/_trimmed.fq.gz/')

## Trimgalore Adjust parameters depending on
## kit used (add clipping and change length for smarter stranded)
if [ -f ${TRIMDIR}/${R1} ]
then
    echo skipping trim for $ID
else
    trim_galore -o $TRIMDIR \
		--length 35 \
		--phred33 \
		--cores 8\
		${FASTQDIR}/${1}
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
	   -U ${TRIMDIR}/${R1}\
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
	     -fq ${TRIMDIR}/${ID}_ribotrim_R1.fq
fi

## Get Trimmed Filename
R1=${ID}_ribotrim_R1.fq

## Run Posttrim FastQC
F1=$(echo $R1 | sed 's/\.fq/_fastqc.html/')
if [ -f ${POSTTRIM_QC}/${F1} ]
then
    echo skipping posttrim QC $R1
else
    fastqc ${TRIMDIR}/${R1} -o ${POSTTRIM_QC}
fi

## Align Trimmed Reads Using Hisat2
if [ -f ${ALIGNDIR}/${ID}_sorted_rf_alignment.bam ]
then
    echo skipping alignment for sample $ID
else
    echo aligning to genome
    echo read1 file ${ID}_ribotrim_R1.fq
    
    hisat2 -p8\
	   --verbose\
	   --phred33\
	   --dta\
	   --fr\
	   --summary-file ${ALIGNDIR}/${ID}_AlignStat.txt\
	   -x $HISAT2_PREFIX\
	   -U ${TRIMDIR}/${ID}_ribotrim_R1.fq\
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
if [ -f ${COUNTDIR}/${ID}_rf_GeneCount.txt ]
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
	     -l 100 -s 15\
	     ${TRIMDIR}/${ID}_ribotrim_R1.fq
fi

# Estimate Transcript Integrity using iteratively generated script
tin.py -i ${ALIGNDIR}/${ID}_sorted_rf_alignment.bam -r $BEDPATH
