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

### Compose "Single Exon BED file for Inner Distance Calculations
# The values in this segment must be determined from the Rn45s nucleotide sequence
# the ribobed variable stores a single "BED" line that defines the entire
# Rn45s nucleotide sequence as one big exon. 
# ribobed=NR_046233.2$"\t"0"\t"13400"\t"NR_046233.2
# ribobed=$ribobed"\t"0"\t"."\t"0"\t"13400"\t"0
# ribobed=$ribobed"\t"1"\t"13400,"\t"0,

echo bed line used for inner distance calculations:
echo -e $ribobed
echo -e $ribobed > ribo.bed

## Estimate Inner distances from each ribosomal alignment
inner_distance.py\
    -i ${ALIGNDIR}/${ID}_sorted_ribo.bam\
    -r ribo.bed\
    -o ${RSEQCDIR}/Ribo_${ID}
rm ribo.bed

## Estimate GC Content of Ribosomal Alignments
read_GC.py\
    -i ${ALIGNDIR}/${ID}_sorted_ribo.bam\
    -o ${RSEQCDIR}/Ribo_${ID}

## Calculate Summary Alignment statistics


## Estimate Fraction of Ribosomal Alignments in Genomic Alignments
genomic=${ID}_sorted_alignment.bam
ribosomal=${ID}_sorted_ribo.bam

gawk -v FS="\t"\
     'BEGIN {rc=0; fc=0}
     (rc > FNR)\
     {
        fc ++
        if(fc == 1)
        {
           unl[$1]
        }
        else
        {
           rib[$1]
        }
     }
    {rc = FNR - 1}
    (fc == 0){aln[$1]}	
    (fc == 1){unl[$1]}	
    (fc == 2){rib[$1]}	
    
    END\
    {
        alnCount = 0
        unlCount = 0
        ribCount = 0
        nulCount = 0

        for(read in rib)
        {
            ribCount ++
            if(read in aln)
            {
               alnCount ++
            }
            else if(read in unl)
            {
               unlCount ++
            }
            else
            {
               nulCount ++
            }
        }
        print ribCount "\t" alnCount "\t" unlCount "\t" nulCount "\t"
    }'\
     <(samtools view -f67 ${ALIGNDIR}/${genomic})\
     <(samtools view -f4 ${ALIGNDIR}/${genomic})\
     <(samtools view -f67 -F256 ${ALIGNDIR}/${ribosomal})\
     >${RSEQCDIR}/${ID}_ribosomal_reads_in_genomic.txt
