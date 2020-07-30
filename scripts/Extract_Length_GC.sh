#!/bin/bash

export SPLITDIR="$(pwd)/Separate_Genes"
export igtf="gencode.vM25.basic.annotation.gtf"
export fasta="GRCm38.primary_assembly.genome.chr.fa"
export SPLITERR=${SPLITDIR}/bad_genes.txt
export RESFILE=Gencode_M25_Basic_Length_GC.tsv

if [ ! -d ${SPLITDIR} ]
then
    echo $SPLITDIR
    mkdir $SPLITDIR
else
    echo $SPLITDIR
    rm -rf $SPLITDIR
    mkdir $SPLITDIR
fi

# Count the number of unique gene entries
genes=$(cat ${igtf} | cut -f9 \
	    | cut -d ";" -f1 \
	    | sort \
	    | uniq \
	    | sed 's/gene_id[[:space:]]//g' \
	    | sed 's/\"//g')
ngenes=$(echo ${genes[@]} | wc -l | tr -d '[:space:]')

echo Unique Gene ID\'s ${ngenes}

# -v SD=${SPLITDIR} Assign output directory path to awk variable 'SD'
# -v FS="\t| |;" Define field delimiters as tab, or single space or semicolon

gawk -v SD=${SPLITDIR} \
     -v FS="\t| |;" \
'{
    gsub(/;\040+/,";");          # Convert "; " delimiters to ";"
    gsub(/\042/, "", $10);       # Strip " quotes from field 10
    arr[NR] = $10;               # Store Field 10 values in an array         
}

# Write exons from input gtf file to separate gtf files named by gene ID
{
    $10 = "\042" $10 "\042";                                # Restore $10 quotes
    if($3 == "exon"){
	print $1"\t"$2"\t"$3"\t"$4"\t"$5		\
	    $6"\t"$7"\t"$8"\t" arr[NR] >> SD"/"arr[NR]".gtf"
    }
}

{
    gsub(";","; ")   
}' $igtf

# Iterate over gene level files, sort by chromosome and start, and merge exons
for f in $(ls ${SPLITDIR}); do
    g=$(echo $f | sed 's/\.gtf//')
    nc=$(cat ${SPLITDIR}/$f | cut -f1 | sort | uniq | wc -l| sed 's/[[:space:]]//g')
    ns=$(cat ${SPLITDIR}/$f | cut -f6 | sort | uniq | wc -l| sed 's/[[:space:]]//g')
    if [ $nc -gt 1 ]; then
	echo $g has features on multiple chromosomes >> ${SPLITERR}
	cat ${SPLITDIR}/$f >> ${SPLITERR}
	continue
    fi

    if [ $ns -gt 1 ]; then
	echo $g has features on multiple strands >> ${SPLITERR}
	cat ${SPLITDIR}/$f >> ${SPLITERR}
	continue
    else
	cat ${SPLITDIR}/$f \
        | sort -k4,4n \
	| bedtools merge -i stdin -c 6,8 -o distinct \
	| gawk -v OFS="\t" '{print $1, $2, $3, $5, ".", $4, $3 - $2}'# \
	| bedtools nuc -fi ${fasta} -bed stdin \
	| gawk '(NR > 1) {print $0}' \
	| bedtools groupby -g 4 -c 10,11,12,13,14,15,16 -o sum \
	| gawk -v OFS="\t" '{print $1, $8, ($3 + $4) / $8}'  # >> ${RESFILE}
    fi
done
