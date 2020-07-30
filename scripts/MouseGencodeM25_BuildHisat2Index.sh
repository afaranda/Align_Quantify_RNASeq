#!/bin/bash
#SBATCH --job-name=BuildM25
#SBATCH --ntasks=8
#SBATCH --mem=512000

## Fetch Genome And GTF
URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25"
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
GTF="gencode.vM25.basic.annotation.gtf.gz"
GEN="GRCm38.primary_assembly.genome.fa.gz"

curl ${URL}/${GTF} --output $GTF
gunzip $GTF
GTF="gencode.vM25.basic.annotation.gtf"

echo ${URL}/${GEN}
curl ${URL}/{$GEN} --output $GEN
gunzip $GEN
GEN="GRCm38.primary_assembly.genome.fa"

## Get Rid of Patches and Scaffolds
BAS=$(echo $GEN | sed 's/genome/genome.chr/')
gawk -v RS=">" -v ORS=">" '(NR < 24) {print $0}' $GEN | head -n -1 > $BAS


## Check Hisat2 Installation
HISAT2_BUILD_EXE=./hisat2-build
if [ ! -x "$HISAT2_BUILD_EXE" ] ; then
        if ! which hisat2-build ; then
                echo "Could not find hisat2-build in current directory or in PATH"
                exit 1
        else
                HISAT2_BUILD_EXE=`which hisat2-build`
	fi
fi

HISAT2_SS_SCRIPT=./hisat2_extract_splice_sites.py
if [ ! -x "$HISAT2_SS_SCRIPT" ] ; then
        if ! which hisat2_extract_splice_sites.py ; then
                echo "Could not find hisat2_extract_splice_sites.py in current directory or in PATH"
                exit 1
        else
                HISAT2_SS_SCRIPT=`which hisat2_extract_splice_sites.py`
        fi
fi

HISAT2_EXON_SCRIPT=./hisat2_extract_exons.py
if [ ! -x "$HISAT2_EXON_SCRIPT" ] ; then
        if ! which hisat2_extract_exons.py ; then
                echo "Could not find hisat2_extract_exons.py in current directory or in PATH"
		exit 1
        else
                HISAT2_EXON_SCRIPT=`which hisat2_extract_exons.py`
        fi
fi

## Extract Exons and Splice Sites
if [ ! -f genome.ss ] ; then
       ${HISAT2_SS_SCRIPT} ${GTF} > genome.ss
       ${HISAT2_EXON_SCRIPT} ${GTF} > genome.exon
fi

## Build Genome
echo hisat2-build -f -p 8 --ss genome.ss --exon genome.exon $BAS M25
hisat2-build -f -p 8 --ss genome.ss --exon genome.exon $BAS M25

hisat2-inspect -s M25

