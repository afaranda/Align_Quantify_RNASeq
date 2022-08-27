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

    
    # eid=1;
    
    # while($eid != "exon_id" && eid < NF){
    # 	eid++;
    # }
    
    # ## Add Exon Row to array
    # gid++;
    # eid++;

    # ## Index exon records by exon ID
    # exon_row[$eid]=$0;

    # ## Index exons on a chromosome
    # chrom[$1];
    # exon_chrom[$eid]=$1;

    # ## Index exons in a gene
    # gene[$gid];
    # exon_gene[$eid]=$gid;

    # ## Index exons in a gene by start position
    # exon_start[$eid]=$4;
    
}

END {
   ## Iterate over chromosomes
   ##   Iterate over gene id's
   ##       Write rows sorted by start to temp file

   ## Call bedtools merge on the temp file and
   ## append output to result file

   ##
#    for(chrom in exons){
#        for(gene in exons[chrom]){
# 	   for(eid in exons[chrom][gene]){
# #	       X=exons[chrom][gene];
# 	       for(eid in exons[chrom][gene]){
# 		   X[0]=0
# 		   X=exons[chrom][gene];
# 		   print(chrom, gene, eid, X[eid]);
# 	       }
# 	   }
#        }
#    }






}
