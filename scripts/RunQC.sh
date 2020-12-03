#!/bin/bash
#SBATCH --job-name=RunQC
#SBATCH --ntasks=4
#SBATCH --mem=32000

#Move slurm output to Alignment Directory
JOBS=$(echo $JOBS | sed 's/afterok://g')

IFS="," read -ra jobs <<< "$JOBS"
for j in ${jobs[@]}; do
    ID=$(head -n 1 slurm-${j}.out)
    mv slurm-${j}.out ${ALIGNDIR}/${ID}_AlignStat.txt
done

ROBS=$(echo $ROBS | sed 's/afterok://g')
IFS="," read -ra robs <<< "$ROBS"
for j in ${robs[@]}; do
    ID=$(head -n 1 slurm-${j}.out)
    mv slurm-${j}.out ${ALIGNDIR}/${ID}_RiboAlignStat.txt
done

# Run Rseqc modules
$(pwd)/scripts/rseqc_modules.sh

# Aggregate QC Results for Main Analysis
multiqc -n $(pwd)/${PWD##*/}_multiqc\
	--config $(pwd)/scripts/multiqc_config.yaml .

# Aggregate QC Results for Ribosomal Analysis
multiqc -n $(pwd)/${PWD##*/}_ribo_multiqc\
	--config $(pwd)/scripts/multiqc_ribo.yaml .


# Aggregate Ribosomal Content estimates
OUTFILE=${RSEQCDIR}/ribosomal_alignment_fractions.txt
echo -e "sample\ttotal\tonly_genomic\tgenomic_and_ribo\tonly_ribo\tno_alignment\tribo_check\tgenomic_check" >> $OUTFILE
for f in $(find ${RSEQCDIR} -regex ".*ribosomal_reads_in_genomic.*")
do
    cat $f >> $OUTFILE
done
