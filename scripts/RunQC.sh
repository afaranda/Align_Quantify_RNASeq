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
multiqc -n $(pwd)/${PWD##*/}_multiqc .
