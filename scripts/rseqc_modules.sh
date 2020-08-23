#!/bin/bash
#SBATCH --job-name=RSeQC_Modules
#SBATCH --mem=32000
#SBATCH --ntasks=4

# Run RSeQC Inner Distance and Gene Body Coverage for Pax6 cKO
# and Wildtype Lens Fiber Cells

export ALIGNDIR=$(pwd)/Alignments
export RSEQCDIR=$(pwd)/RSeQC_Results
export ALLBEDPATH=/work/abf/MouseEnsembl100/rseqc_gene_models.bed
export ANALYSISID=Ens100_Fibers_All_Genes

# Create Directory for RSeQC results if it doesn't exist
if [ ! -d $RSEQCDIR ]
then
    mkdir $RSEQCDIR
fi

# Generate a comma separated list of all target bam files
bf=""
echo Processing BAM Files:
for b in $(find $ALIGNDIR | grep LF | grep bam$)
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
	-r ${ALLBEDPATH}\
	-o ${RSEQCDIR}/${ANALYSISID}_${fn}
done

# Estimate Read Distributions
for b in $(echo $bf | sed 's/,/ /g')
do
    fn=$(echo $b | sed "s|$ALIGNDIR/||g"| sed 's/_sorted_alignment\.bam//g')
    read_distribution.py\
	-i $b\
	-r ${ALLBEDPATH} > ${RSEQCDIR}/${ANALYSISID}_${fn}
done

# Estimate GC Content for Aligned Reads.
for b in $(echo $bf | sed 's/,/ /g')
do
    fn=$(echo $b | sed "s|$ALIGNDIR/||g"| sed 's/_sorted_alignment\.bam//g')
    read_GC.py\
	-i $b\
	-o ${RSEQCDIR}/${ANALYSISID}_${fn}
done


# Generate BAM file alignment statistics
for b in $(echo $bf | sed 's/,/ /g')
do
    fn=$(echo $b | sed "s|$ALIGNDIR/||g"| sed 's/_sorted_alignment\.bam//g')
    echo $fn > ${RSEQCDIR}/${ANALYSISID}_${fn}_bam_stats.txt
    bam_stat.py\
	-i $b >> ${RSEQCDIR}/${ANALYSISID}_${fn}_bam_stats.txt
done


# Estimate Transcript Integrity
tin.py\
    -i ${ALIGNDIR}\
    -r ${ALLBEDPATH}

# Iterate over tin.py results and move to results RSeQC Results Directory
for f in $(find . -maxdepth 1 -regex .*sorted_alignment.summary.txt)
do
    fn=$(echo $f | sed 's/sorted_alignment/tin/')
    mv $f $(pwd)/${RSEQCDIR}/${fn}
done

for f in $(find . -maxdepth 1 -regex .*sorted_alignment.tin.xls)
do
    fn=$(echo $f | sed 's/_sorted_alignment//')
    mv $f $(pwd)/${RSEQCDIR}/${fn}
done

# Estimate Gene Body Coverage using the specified bed file
geneBody_coverage.py\
    -i ${ALIGNDIR}\
    -r ${ALLBEDPATH}\
    -o ${RSEQCDIR}/${ANALYSISID}
