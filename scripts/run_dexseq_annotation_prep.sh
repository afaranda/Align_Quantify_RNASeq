#!/bin/bash
#SBATCH --mem=32000

python3 dexseq_prepare_annotation.py Mus_musculus.GRCm39.104.gtf EnsMm_104_DEXSeq_Annot.gff
