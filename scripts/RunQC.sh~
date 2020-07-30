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

# 
multiqc .
