#!/bin/bash
#SBATCH --job-name=RSeQC_Modules
#SBATCH --mem=32000
#SBATCH --ntasks=4

#
# TO_DO: Either parallelize the Gene_Body_Coverage Analysis or
#        Figure out a way to run it over a smaller subset of
#        genes -- the whole genome takes way too long

# Uncomment and edit the following variables if running standalone
# export ALIGNDIR=$(pwd)/Alignments
# export RSEQCDIR=$(pwd)/RSeQC_Results
# export BEDPATH=/work/abf/MouseEnsembl100/rseqc_gene_models.bed
export ANALYSISID=Main_All_Genes
export FILTER_BED=0 # Set to 1 if bam files should be filtered

# Create Directory for RSeQC results if it doesn't exist
if [ ! -d $RSEQCDIR ]
then
    mkdir $RSEQCDIR
fi

# Generate a comma separated list of all target bam files
bf=""
echo Processing BAM Files:
for b in $(find $ALIGNDIR -regex .*_sorted_alignment.bam$)
do
    echo $b
    bf="${bf},${b}"
done
#echo $bf

# Estimate Inner Distances from each alignment
for b in $(echo $bf | sed 's/,/ /g')
do
    fn=$(echo $b | sed "s|$ALIGNDIR/||g"| sed 's/_sorted_alignment\.bam//g')
    inner_distance.py\
	-i $b\
	-r ${BEDPATH}\
	-o ${RSEQCDIR}/${ANALYSISID}_${fn}
done

# Estimate Read Distributions
for b in $(echo $bf | sed 's/,/ /g')
do
    fn=$(echo $b | sed "s|$ALIGNDIR/||g"| sed 's/_sorted_alignment\.bam//g')
    read_distribution.py\
	-i $b\
	-r ${BEDPATH} > ${RSEQCDIR}/${ANALYSISID}_${fn}_read_dists.txt
done

# Estimate GC Content for Aligned Reads. If FILTER_BED == 1, then filter the
# BAM file using samtools view to pass only alignments that intersect with
# the bed file specified in BEDPATH.
if [ $FILTER_BED == 1 ]
then
    for b in $(echo $bf | sed 's/,/ /g')
    do
	fn=$(echo $b | sed "s|$ALIGNDIR/||g"| sed 's/_sorted_alignment\.bam//g')
	read_GC.py\
	    -i <(samtools view -b -L ${BEDPATH} $b)\
	    -o ${RSEQCDIR}/${ANALYSISID}_${fn}
		 
    done
else
    for b in $(echo $bf | sed 's/,/ /g')
    do
	fn=$(echo $b | sed "s|$ALIGNDIR/||g"| sed 's/_sorted_alignment\.bam//g')
	read_GC.py\
	    -i $b\
	    -o ${RSEQCDIR}/${ANALYSISID}_${fn}
		 
    done    
fi

# Generate BAM file alignment statistics
for b in $(echo $bf | sed 's/,/ /g')
do
    fn=$(echo $b | sed "s|$ALIGNDIR/||g"| sed 's/_sorted_alignment\.bam//g')
    echo $fn > ${RSEQCDIR}/${ANALYSISID}_${fn}_bam_stats.txt
    bam_stat.py\
	-i $b >> ${RSEQCDIR}/${ANALYSISID}_${fn}_bam_stats.txt
done

# Calculate Fragment Sizes
for b in $(echo $bf | sed 's/,/ /g')
do
    fn=$(echo $b | sed "s|$ALIGNDIR/||g"| sed 's/_sorted_alignment\.bam//g')
    RNA_fragment_size.py\
        -i $b\
	-r ${BEDPATH} >> ${RSEQCDIR}/${ANALYSISID}_${fn}_frag_sizes.txt
done

# Estimate Gene Body Coverage using the specified bed file
geneBody_coverage.py\
    -i ${bf}\
    -r /work/abf/MouseEnsembl101/rseqc_perinatal_lethal_models.bed\
    -o ${RSEQCDIR}/${ANALYSISID}_lethal

geneBody_coverage.py\
    -i ${bf}\
    -r /work/abf/MouseEnsembl101/rseqc_crystallin_transcript_models.bed\
    -o ${RSEQCDIR}/${ANALYSISID}_crystallin
