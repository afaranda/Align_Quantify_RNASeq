#!/bin/bash
#SBATCH --job-name=SE_ribo
#SBATCH --mem=128000
#SBATCH --ntasks=8


## Get Trimmed Filenames and sample label
ID=$(echo $1 | sed 's/_L[0-9]\{3\}_R[12]_[0-9]\{3\}\.fastq\.gz//')
R1=$(echo $1 | sed 's/\.fastq\.gz/_val_1.fq.gz/')

echo $ID

## Align Trimmed Reads To Ribosomal Index Using Hisat2
# if [ -f ${ALIGNDIR}/${ID}_sorted_ribo.bam ]
# then
#    echo skipping alignment for sample $ID
# else   
#     hisat2 -p8\
# 	   --verbose\
# 	   --phred33\
# 	   --dta\
# 	   --fr\
# 	   -x Rn45s_Index\
# 	   -U ${TRIMDIR}/${R1}\
# 	   -S ${ALIGNDIR}/${ID}_ribo_reads.sam

#     ## compress, sort, and index alignments
#     samtools view -@ 8 -bS ${ALIGNDIR}/${ID}_ribo_reads.sam > ${ALIGNDIR}/${ID}_ribo_reads.bam
    
#     samtools sort -@ 8 -m 12G -o ${ALIGNDIR}/${ID}_sorted_ribo.bam ${ALIGNDIR}/${ID}_ribo_reads.bam
    
#     samtools index ${ALIGNDIR}/${ID}_sorted_ribo.bam
    
#     rm ${ALIGNDIR}/${ID}_ribo_reads.bam ${ALIGNDIR}/${ID}_ribo_reads.sam
# fi

## Estimate GC Content of Ribosomal Alignments
read_GC.py\
    -i ${ALIGNDIR}/${ID}_byname_ribo.bam\
    -o ${RSEQCDIR}/Ribo_${ID}

## Estimate Fraction of Ribosomal Alignments in Genomic Alignments
# genomic=${ID}_sorted_alignment.bam
# ribosomal=${ID}_sorted_ribo.bam

# gawk -v FS="\t"\
#      'BEGIN {rc=0; fc=0}
#      (rc > FNR)\
#      {
#         fc ++
#         if(fc == 1)
#         {
#            unl[$1]
#         }
#         else
#         {
#            rib[$1]
#         }
#      }
#     {rc = FNR - 1}
#     (fc == 0){all[$1]}	
#     (fc == 1){gen[$1]}	
#     (fc == 2){rib[$1]}	
    
#     END\
#     {
#         allCount = 0
#         genCount = 0
#         ribCount = 0
#         nulCount = 0
#         inxCount = 0

#         for(read in all)
#         {
#             allCount ++
#             if((read in rib) && (read in gen))
#             {
#                inxCount ++
#             }
#             else if((read in rib) && !(read in gen))
#             {
#                ribCount ++
#             }
#             else if(!(read in rib) && (read in gen))
#             {
#                genCount ++
#             }
#             else
#             {
#                nulCount ++
#             }
#         }
#         ribCheck = 0
#         for(read in rib)
#         {
#             ribCheck ++
#         }
        
#         genCheck = 0
#         for(read in gen)
#         {
#             genCheck ++
#         }
         
#         print allCount "\t" genCount "\t" inxCount "\t" ribCount "\t"\
#               nulCount "\t" ribCheck "\t" genCheck
#     }'\
#      <(samtools view ${ALIGNDIR}/${genomic})\
#      <(samtools view -F4 ${ALIGNDIR}/${genomic})\
#      <(samtools view -F4 ${ALIGNDIR}/${ribosomal})\
#      >${RSEQCDIR}/${ID}_ribosomal_reads_in_genomic.txt

# # genomic=${ID}_sorted_alignment.bam
# # ribosomal=${ID}_sorted_ribo.bam

# # gawk -v FS="\t"\
# #      'BEGIN {rc=0; fc=0}
# #      (rc > FNR)\
# #      {
# #         fc ++
# #         if(fc == 1)
# #         {
# #            unl[$1]
# #         }
# #         else
# #         {
# #            rib[$1]
# #         }
# #      }
# #     {rc = FNR - 1}
# #     (fc == 0){aln[$1]}	
# #     (fc == 1){unl[$1]}	
# #     (fc == 2){rib[$1]}	
    
# #     END\
# #     {
# #         alnCount = 0
# #         unlCount = 0
# #         ribCount = 0
# #         nulCount = 0

# #         for(read in rib)
# #         {
# #             ribCount ++
# #             if(read in aln)
# #             {
# #                alnCount ++
# #             }
# #             else if(read in unl)
# #             {
# #                unlCount ++
# #             }
# #             else
# #             {
# #                nulCount ++
# #             }
# #         }
# #         print ribCount "\t" alnCount "\t" unlCount "\t" nulCount "\t"
# #     }'\
# #      <(samtools view -f67 ${ALIGNDIR}/${genomic})\
# #      <(samtools view -f4 ${ALIGNDIR}/${genomic})\
# #      <(samtools view -f67 -F256 ${ALIGNDIR}/${ribosomal})\
# #      >${RSEQCDIR}/${ID}_ribosomal_reads_in_genomic.txt
# # Subtract Ribosomal Reads from Genomic Alignments

# echo "removing reads aligned in" ${ribosomal} "from" ${genomic}
# echo 

# gawk '\
#      (NR == FNR) {tag[$1]; next}
#      (NR != FNR && $1 ~/^@/){print $0; next}
#      (NR != FNR && !($1 in tag)){print $0}
#      '\
#      <(samtools view -F4 ${ALIGNDIR}/${ribosomal})\
#      <(samtools view -h ${ALIGNDIR}/${genomic})\
#      | samtools view -b -o ${ALIGNDIR}/${ID}_ribo_filtered_alignment.bam

# # Index filtered bam File
# samtools index ${ALIGNDIR}/${ID}_ribo_filtered_alignment.bam

# # Generate Counts from Ribo-stripped alignments
# if [ -f ${COUNTDIR}/${ID}_RiboFilterGeneCount.txt ]
# then
#     echo skipping htseq-count for ${ID}_ribo_filtered_alignment.bam
# else
#     echo counting ${ID}_ribo_filtered_alignment.bam
#     htseq-count \
# 	-i gene_id -r pos -f bam -s reverse -m union --type exon \
# 	${ALIGNDIR}/${ID}_ribo_filtered_alignment.bam \
# 	$GTFPATH > ${COUNTDIR}/${ID}_RiboFilterGeneCount.txt
# fi



