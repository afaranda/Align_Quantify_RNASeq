#!/bin/bash
###############################################################################
#
# File:    Gencode_Extract_EU_PL_Length_GC.sh
#
# Purpose: Calculate gene level estimates of length and GC content using
#          using features defined in a Gencode formatted GTF File
#
#          This script implements two methods for length and gc content
#          estimation.  The exon-union method (EU) and the principal-longest
#          (PL) method.
#
#          The EU method assumes all annotated exons contribute to the
#          transcript length. In the EU method, overlapping exons are merged
#          into a "union-exon, the sum of all exon lengths is considered the
#          gene level length, and the sum of all G and C nucleotides in each
#          exon is used to estimate a gene level gc content.
#
#          The PL method selects a single, representative transcript for any
#          gene annotated with more than one transcript. If there is one or
#          more transcripts bearing the "appris_principal_1" tag, the longest
#          of these will be selected to represent the gene.  If there are no
#          transcripts bearing this tag, then the longest overall transcript
#          is selected to represent the gene.
#
#          In the Gencode GTF attribute string (field 9 on tab delimiter)
#          gene_id is allways the first attribute, and for exons and
#          transcripts, transcript_id is allways the second attribute.
#
#         
#
# Author:  Adam Faranda
# Created: October 5, 2020
#
###############################################################################

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
        i=0
	printf $1"\t"$2"\t"$3"\t"$4"\t"$5\
	    "\t"$6"\t"$7"\t"$8"\t" >> SD"/"arr[NR]".gtf"
        for(i=9; i <=NF; i++){
            if(i % 2 == 1){
                printf $i" " >> SD"/"arr[NR]".gtf"
            }
            else{
                printf $i"; " >> SD"/"arr[NR]".gtf"
            }
        }
        print "" >> SD"/"arr[NR]".gtf"
    }
}

{
    gsub(";","; ")   
}' $igtf

# Iterate over gene level files, calculate exon union length and GC content,
# calculate primary transcript length and GC content according to Elkon's
# rules, tabulate gene level length / GC statistics
for f in $(ls ${SPLITDIR}); do
    
    g=$(echo $f | sed 's/\.gtf//')        # Store the Gene ID in a variable

    # Check the number of chromsomes (should only be one)
    nc=$(
	cat ${SPLITDIR}/$f\
	    | cut -f1 | sort | uniq\
	    | wc -l| sed 's/[[:space:]]//g'
      )

    # Check the number of strand orientations (should only be one)
    ns=$(
	cat ${SPLITDIR}/$f\
	    | cut -f6 | sort | uniq\
	    | wc -l| sed 's/[[:space:]]//g'
      )

    # Stop processing this gene on chromsome error
    if [ $nc -gt 1 ]; then
	echo $g has features on multiple chromosomes >> ${SPLITERR}
	cat ${SPLITDIR}/$f >> ${SPLITERR}
	continue
    fi

    # Stop processign this gene on strand error, otherwise process
    if [ $ns -gt 1 ]; then
	echo $g has features on multiple strands >> ${SPLITERR}
	cat ${SPLITDIR}/$f >> ${SPLITERR}
	continue

    else
	# read set of exons, sort by start position, calculate exon lengths
	# and nucleotide content -- passed to gawk for processing
	cat ${SPLITDIR}/$f \
            | sort -k4,4n\
	    | bedtools nuc -s -fi ${fasta} -bed stdin\
	    | gawk \
	    '(NR > 1 && NR == FNR)\
             {
                 FS="\t"                                 # split on tab
                 gsub(/;\040+/,";")                      # "; " becomes ";"
                 gsub(/\042/,"")                         # drop quotes
                 split($9,atr,";| ")               # Store field 9 in array atr      
                 gid[atr[2]]                       # Store gene id as array key
                 tlen[atr[4]]=tlen[atr[4]] + $18   # sum exon lengths by txid
                 if($0 ~ /appris_principal_1/){    # Flag "P" if principal
                    tprn[atr[4]]="P"                 
                 }
                 else{
                     tprn[atr[4]]="N"              # Flag "N" if not principal
                 }
                 txgc[atr[4]]=txgc[atr[4]] + $13 + $14
             }
             # After loading 
             (FNR != NR && $1 in gid){
                    FS="\t"                        
                    for(g in gid){
                        l = 0                      # max non-principal tx len
                        m = 0                      # max principal tx len 
                        longest = ""               # id of longest "N" tx  
                        principal =""              # id of longest "P" tx
                        
                        # Iterate over transcripts and identify longest
                        # principal and non-principal
                        for(t in tlen){
                            if(tprn[t] == "P" && tlen[t] > m){  
                                principal=t
                                m=tlen[t]
                            }
                            else if(tlen[t]>l && m == 0){
                                longest=t
                                l=tlen[t]
                            }
                            
                        }

                        # Return either the longest principal transcript
                        if(m > 0)
                            print $1"\t" principal "\t" $2 "\t" $3 "\t"\
                                  tlen[principal]"\t"\
                                  txgc[principal] / tlen[principal]"\t"\
                                  tprn[principal]
                        
                        # Or if none exists, the longest transcript 
                        else
                            print $1"\t" longest "\t" $2 "\t" $3 "\t"\
                                  tlen[longest]"\t"\
                                  txgc[longest] / tlen[longest] "\t"\
                                  tprn[longest]
                    }
                }' /dev/stdin\
		    <(
		# Calculate Exon-Union Length & GC content using
		# bedtools merge to flatten overlapping exons
		       cat ${SPLITDIR}/$f\
			   | sort -k4,4n\
			   | gawk -v FS="\t| |;"\
           	           '\
                           {
                               gsub(/;\040+/,";")
                               gsub(/\042/,"")
                               for(i=1; i<9;i++){printf($i"\t")}
                               print($10)
                           }'\
			   | bedtools merge -s -i stdin -c 7,9 -o distinct\
			   | gawk -v OFS="\t"\
				  '{print $1, $2, $3, $5, ".", $4, $3 - $2}'\
			   | bedtools nuc -s -fi ${fasta} -bed stdin\
			   | gawk '(NR > 1)'\
		       	   | bedtools groupby -g 4\
				      -c 10,11,12,13,14,15,16 -o sum\
			   | gawk -v OFS="\t" '{print $1, $8, ($3 + $4) / $8}'
	            ) >> ${RESFILE}
    fi
done
