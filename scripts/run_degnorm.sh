#!/bin/bash
#SBATCH --job-name=degnorm
#SBATCH --mem=128000
#SBATCH --ntasks=8

# Call R script encapsulating DegNorm process
Rscript scripts/DegNorm_Estimate_Degradation.R $1

