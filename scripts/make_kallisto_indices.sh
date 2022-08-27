#!/bin/bash
#SBATCH --job-name=BUILD_MOUSE_INDEX
#SBATCH --ntasks=8                                                                                                                                                
#SBATCH --mem=128000

BEDFILE=rseqc_gene_models.bed
GENOME=Mus_musculus.GRCm39.dna.primary_assembly.fa
INDEX_ID=EnsMm_104_Kallisto

# Build Kallisto Index for Exons using bed file
# generated by ensembl_make_rseqc_bedfile.sh

kallisto index -i ${INDEX_ID}_exon\
	 <(\
	   bedtools getfasta\
		    -split\
		    -nameOnly\
		    -fi $GENOME\
                    -bed $BEDFILE\
    )

# Build Kallisto Index for total transcript using bed file
# generated by ensembl_make_rseqc_bedfile.sh

kallisto index -i ${INDEX_ID}_total\
	 <(\
	   bedtools getfasta\
		    -nameOnly\
		    -fi $GENOME\
                    -bed $BEDFILE\
    )