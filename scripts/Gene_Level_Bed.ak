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

    ## Iterate over fields to find transcript_id
    tid=1;
    while ($tid != "transcript_id"){
	tid++;
    }
    tid++;

    ## Iterate over fields to find exon_id
    eid=1;
    while ($eid != "exon_id"){
	eid++;
    }
    eid++;

    gsub("\042","",$gid)
    gsub("\042","",$tid)
    gsub("\042","",$eid)
    # Add exon, indexed by gene and record number
    exons_g[$gid][NR]=$1"\t"$4"\t"$5"\t"$gid"\t"".""\t"$7;

    # Index start and end positions of exons      
    ex_start[$tid][$eid]=$4;
    ex_end[$tid][$eid]=$5;
    
}

($3 ~ /transcript/){
    ## Iterate over fields to find transcript_id
    tid=1;
    while ($tid != "transcript_id"){
	tid++;
    }
    tid++;

    gsub("\042","",$tid)
    ## Store transcript record
    tx[$tid]=$1"\t"$4 - 1"\t"$5"\t"$tid"\t"".""\t"$7;
}

($3 ~ /gene/){
    ## Iterate over fields to find gene_id
    gid=1;
    while($gid != "gene_id"){
	gid++;	  
    }
    gid++;
    
    gsub("\042","",$gid)
    ## Store gene record
    gn[$gid]=$1"\t"$4"\t"$5"\t"$gid"\t"".""\t"$7;
}


END {

    ## Write collapsed gene-level BED file to store flattened exons
    for(g in exons_g){
        for(r in exons_g[g]){
	    printf exons_g[g][r]"\n" >> "temp.bed";
	}
	system("bedtools sort -i temp.bed > sorted.bed");
	call="bedtools merge -s -o distinct -c 4,5,6 -i sorted.bed >> result.bed"
	system(call);
	close("temp.bed");
	close("sorted.bed");
	system("rm temp.bed");
	system("rm sorted.bed");
    }
    close("result.bed");
    
    ## Write Transcripts to a file in BED12 format
    for(tid in ex_start){
     
	# Sort exons of a transcript in ascending numerical order 
	m=asort(ex_start[tid],exs,"@val_num_asc");
	n=asort(ex_end[tid],exe,"@val_num_asc");
	
	# Split fields from the corresponding transcript record
	o=split(tx[tid],tx_row,"\t");
	
	# Write fields to BED12 file
	for(i=1; i < 4; i++){
	    printf tx_row[i]"\t" >> "rseqc_transcripts.bed";
	}
	printf tid"\t""0""\t"tx_row[6]"\t" >> "rseqc_transcripts.bed";
	
	for(i=2; i < 4; i++){
	    printf tx_row[i]"\t" >> "rseqc_transcripts.bed";
	}
	printf "0""\t"m"\t" >> "rseqc_transcripts.bed";
	for(i=1; i <=m; i++){
	    printf (exe[i] - (exs[i] - 1))"," >> "rseqc_transcripts.bed";
	}
	printf("\t") >> "rseqc_transcripts.bed";
	for(i=1; i <=m; i++){
	    printf ((exs[i] -1) - tx_row[2])"," >> "rseqc_transcripts.bed";
	}
	printf("\n")  >> "rseqc_transcripts.bed";
    }
    close("rseqc_transcripts.bed");

    ## Load "gene level collapsed exons into an array"
    flag=(getline line < "result.bed");
    idx=0;
    print(line)
    while(flag == 1){
	print("line in loop "line);
	idx++;
	gene_bed[idx] = line;
	flag=(getline line < "result.bed");
    }

    ## Index start and stop positions from gene level bed
    for(i=1; i <= idx; i++){
	split(gene_bed[i],rec,"\t");
	gx_start[rec[4]][i]=rec[2];
	gx_end[rec[4]][i]=rec[3];
    }

    ## Write Collapsed genes to a file in BED12 format
    for(gid in gx_start){
     
	# Sort exons of a transcript in ascending numerical order 
	m=asort(gx_start[gid],gxs,"@val_num_asc");
	n=asort(gx_end[gid],gxe,"@val_num_asc");
	
	# Split fields from the corresponding transcript record
	o=split(gn[gid],gn_row,"\t");
	
	# Write fields to BED12 file
	for(i=1; i < 4; i++){
	    printf gn_row[i]"\t" >> "rseqc_genes.bed";
	}
	printf gid"\t""0""\t"gn_row[6]"\t" >> "rseqc_genes.bed";
	
	for(i=2; i < 4; i++){
	    printf gn_row[i]"\t" >> "rseqc_genes.bed";
	}
	printf "0""\t"m"\t" >> "rseqc_genes.bed";
	for(i=1; i <=m; i++){
	    printf (gxe[i] - gxs[i])"," >> "rseqc_genes.bed";
	}
	printf("\t") >> "rseqc_genes.bed";
	for(i=1; i <=m; i++){
	    printf (gxs[i]  - gn_row[2])"," >> "rseqc_genes.bed";
	}
	printf("\n")  >> "rseqc_genes.bed";
    }
    close("rseqc_genes.bed");
}
