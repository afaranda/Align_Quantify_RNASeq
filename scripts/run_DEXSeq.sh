#!/bin/bash
#SBATCH --mem=64000
R CMD BATCH scripts/DEXSeq_Differential_Exons.R DEXSeq.Rout
