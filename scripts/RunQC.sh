#!/bin/bash
#SBATCH --job-name=RunQC
#SBATCH --ntasks=4
#SBATCH --mem=32000
#RSEQCDIR=$(pwd)/RSeQC_Results
ANALYSISID=Main_All_Genes

#Move slurm output to Alignment Directory
#JOBS=$(echo $JOBS | sed 's/afterok://g')

#IFS="," read -ra jobs <<< "$JOBS"
#for j in ${jobs[@]}; do
#    ID=$(head -n 1 slurm-${j}.out)
#    mv slurm-${j}.out ${ALIGNDIR}/${ID}_AlignStat.txt
#done

#ROBS=$(echo $ROBS | sed 's/afterok://g')
#IFS="," read -ra robs <<< "$ROBS"
#for j in ${robs[@]}; do
#    ID=$(head -n 1 slurm-${j}.out)
#    mv slurm-${j}.out ${ALIGNDIR}/${ID}_RiboAlignStat.txt
#done

# Run Rseqc modules
$(pwd)/scripts/${RSEQCMOD}

# Aggregate QC Results for Main Analysis
multiqc -n $(pwd)/${PWD##*/}_multiqc\
	--config $(pwd)/scripts/multiqc_config.yaml .

# Aggregate QC Results for Ribosomal Analysis
multiqc -n $(pwd)/${PWD##*/}_ribo_multiqc\
	--config $(pwd)/scripts/multiqc_ribo.yaml .

# Iterate over tin.py results and move to results RSeQC Results Directory
for f in $(find . -maxdepth 1 -regex .*sorted_rf_alignment.summary.txt)
do
    fn=$(echo $f | sed 's/sorted_rf_alignment/tin/'| sed 's/^\.\///')
    mv $f ${RSEQCDIR}/${ANALYSISID}_${fn}
done

for f in $(find . -maxdepth 1 -regex .*sorted_rf_alignment.tin.xls)
do
    fn=$(echo $f | sed 's/_sorted_rf_alignment//'| sed 's/^\.\///')
    mv $f ${RSEQCDIR}/${ANALYSISID}_${fn}
done

## Launch DegNorm
sbatch run_degnorm.sh ribo

# # Aggregate Ribosomal Content estimates
# OUTFILE=${RSEQCDIR}/ribosomal_alignment_fractions.txt
# echo -e "sample\ttotal\tonly_genomic\tgenomic_and_ribo\tonly_ribo\tno_alignment\tribo_check\tgenomic_check" >> $OUTFILE
# for f in $(find ${RSEQCDIR} -regex ".*ribosomal_reads_in_genomic.*")
# do
#     cat $f >> $OUTFILE
# done
