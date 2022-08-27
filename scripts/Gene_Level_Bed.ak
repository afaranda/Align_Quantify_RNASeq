#!/usr/bin/gawk -f
###################################################################
#
# File: Gene_Level_Bed.ak
# Purpose: Collapse exons in GTF file from transcripts to genes
#          and store resulting intervals in a BED file
# Created: August 24, 2022
# Author: Adam Faranda
#
###################################################################




## Set delimiters and environment variables
BEGIN {
    FS="\t| |; ";
}

## At each record with the value "exon" in field 3
($3 ~ /exon/){

    ## Populate an array indexed by chromosome, then by gene id and
    ## then by start position to capture GTF file rows corresponding
    ## to each gene
    
    ## Iterate over fields to find gene_id
    gid=1;
    while($gid != "gene_id"){
	gid++;	  
    }
    gid++;

    # Add exon, indexed by gene and record number
    exons[$gid][NR]=$0;    
}


END {
    for(g in exons){
        for(r in exons[g]){
	    printf exons[g][r]"\n" >> "temp.bed";
	}
	system("bedtools sort -i temp.bed > sorted.bed");
	system("bedtools merge -o distinct -c 7 -i sorted.bed >> result.bed");
	close("temp.bed");
	close("sorted.bed");
	system("rm temp.bed");
	system("rm sorted.bed");
    }
}
