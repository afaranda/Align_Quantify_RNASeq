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

# Generate Bed File
if [ -f ${HISAT2_INDEXES}/rseqc_gene_models.bed ]; then
    awk'{if($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"}' ${GTFPATH}\
	| gtf2bed - > rseqc_gene_models.bed
fi 

multiqc .
