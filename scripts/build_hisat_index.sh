#!/bin/sh                                                                                                                                        
#SBATCH --job-name=BUILD_MOUSE_INDEX
#SBATCH --ntasks=8                                                                                                                                                
#SBATCH --mem=512000

# Downloads sequence for the GRCm38 release 96 version of M. musculus (mouse) from                                                               
# Ensembl.                                                                                                                                       
#                                                                                                                                                
# By default, this script builds and index for just the base files,                                                                              
# since alignments to those sequences are the most useful.  To change                                                                            
# which categories are built by this script, edit the CHRS_TO_INDEX                                                                              
# variable below.                                                                                                                                
#                                                                                                                                                

export PATH=$PATH:/home/abf/bin
#declare -a CHR=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y MT)
#declare -a GENOME=()

ENSEMBL_RELEASE=100
ENSEMBL_GRCm38_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/dna
ENSEMBL_GRCm38_GTF_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/mus_musculus
#GTF_FILE=Mus_musculus.GRCm38.${ENSEMBL_RELEASE}.chr.gtf # Excludes unplaced contigs
GTF_FILE=Mus_musculus.GRCm38.${ENSEMBL_RELEASE}.gtf

get() {
        file=$1
	if ! wget --version >/dev/null 2>/dev/null ; then
                if ! curl --version >/dev/null 2>/dev/null ; then
                        echo "Please install wget or curl somewhere in your PATH"
                        exit 1
                fi
                curl -o `basename $1` $1
                return $?
        else
                wget -nv $1
                return $?
	fi
}

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

#rm -f genome.fa
# Un comment this block if retrieving individual chromosomes
#for c in ${CHR[@]}; do
#    
#    F="Mus_musculus.GRCm38.dna.chromosome.$c.fa"
#    G=$(echo $F | sed 's/Mus_musculus\.GRCm38\.dna\.chromosome\./chr./')
#    if [ ! -f $G ] ; then
#        get ${ENSEMBL_GRCm38_BASE}/$F.gz || (echo "Error getting $F" && exit 1)
#        gunzip $F.gz || (echo "Error unzipping $F" && exit 1)
#        mv $F "chr.$c.fa"
#    fi
#    GENOME=("${GENOME[@]}" "chr.$c.fa")
#    
#done

#GENOME=$(echo ${GENOME[@]} | sed 's/\s/,/g')

if [ ! -f $GTF_FILE ] ; then
       get ${ENSEMBL_GRCm38_GTF_BASE}/${GTF_FILE}.gz || (echo "Error getting ${GTF_FILE}" && exit 1)
       gunzip ${GTF_FILE}.gz || (echo "Error unzipping ${GTF_FILE}" && exit 1)
fi

if [ ! -f genome.ss ] ; then
       ${HISAT2_SS_SCRIPT} ${GTF_FILE} > genome.ss
       ${HISAT2_EXON_SCRIPT} ${GTF_FILE} > genome.exon
fi



#Un-comment if retrieving complete genome (with patches)
F="Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
GENOME=$(echo $F | sed 's/\.gz//')

if [ ! -f $GENOME ] ; then
    get ${ENSEMBL_GRCm38_BASE}/${F} || (echo "Error Fetching $F" && exit 1)
    gunzip $F
fi

#CMD="${HISAT2_BUILD_EXE} -p 8 $GENOME --ss genome.ss --exon genome.exon genome_tran"
#echo Running $CMD
#if $CMD ; then
#    echo "genome index built; you may remove fasta files"
#else
#    echo "Index building failed; see error message"
#fi
echo hisat2-build -f -p 8 --ss genome.ss --exon genome.exon $GENOME genome_tran
hisat2-build -f -p 8 --ss genome.ss --exon genome.exon $GENOME genome_tran
Rscript /home/abf/bin/flatten_exon.R
