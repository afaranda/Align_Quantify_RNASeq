#!/bin/bash
#SBATCH --job-name=RSeQC_Modules
#SBATCH --mem=32000
#SBATCH --ntasks=4

export ALIGNDIR=$(pwd)/Alignments
export RSEQCDIR=$(pwd)/RSeQC_Results
export ALLBEDPATH=/work/abf/MouseEnsembl100/rseqc_gene_models.bed
export ANALYSISID=EnsMm100_All

# Create Directory for RSeQC results if it doesn't exist
if [ ! -d $RSEQCDIR ]
then
    mkdir $RSEQCDIR
fi


# Generate a comma separated list of all target bam files
bf=""
for b in $(find $ALIGNDIR | grep bam$)
do
    bf="${bf},${b}"
done
echo $bf

# Estimate Gene Body Coverage using the specified bed file
geneBody_coverage.py\
    -i ${ALIGNDIR}\
    -r ${ALLBEDPATH}\
    -o ${RSEQCDIR}/${ANALYSISID}

# Estimate Inner Distances from each alignment
for b in $(echo $bf | sed 's/,/ /g')
do
    inner_distance.py\
	-i $b\
	-r ${ALLBEDPATH}\
	-o ${RSEQCDIR}/${ANALYSISID}
done
