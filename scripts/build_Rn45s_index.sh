#!/bin/sh                                                                                                                                        
#SBATCH --job-name=BUILD_MOUSE_INDEX
#SBATCH --ntasks=8                                                                                                                                                
#SBATCH --mem=512000

# Build a HISAT2 Index for the Human Ribosomal RNA Pre-Spliced Transcript
# Retrieved from NCBI on May 11, 2021 (File RNA45SN5.fasta)


hisat2-build -f -p 8 RNA45SN5.fasta RNA45SN5_Index

