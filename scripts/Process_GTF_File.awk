#!/usr/bin/gawk -f
###################################################################
#
# File: Process_GTF_File.awk
# Purpose: Collapse exons in GTF file from transcripts to genes
#          and store resulting intervals in a BED file
# Created: August 24, 2022
# Author: Adam Faranda
#
###################################################################




## Set delimiters and environment variables
BEGIN {
    FS="\t| |; ";
    genome="/work/abf/MouseEnsembl104";
    genome=genome"/Mus_musculus.GRCm39.dna.primary_assembly.fa";
}

## At each record with the value "exon" in field 3
($3 ~ /exon/){
    
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

    ## Strip quotes from ensembl ID's
    gsub("\042","",$gid)
    gsub("\042","",$tid)
    gsub("\042","",$eid)
    
    ## Store exons in an array, indexed by gene ($gid) and record number (NR)
    ## Each exon is represented by a tab delimited string
    exons_g[$gid][NR]=$1"\t"$4 - 1"\t"$5"\t"$gid"\t"".""\t"$7;

    ## Index start and end positions of exons      
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

    ## Store transcript record as a tab delimited string
    gsub("\042","",$tid)
    tx[$tid]=$1"\t"$4 - 1"\t"$5"\t"$tid"\t"".""\t"$7;
}

($3 ~ /gene/){
    ## Iterate over fields to find gene_id
    gid=1;
    while($gid != "gene_id"){
	gid++;	  
    }
    gid++;
    
    ## Store gene record as a tab delimited string
    gsub("\042","",$gid)
    gn[$gid]=$1"\t"$4"\t"$5"\t"$gid"\t"".""\t"$7;
}


END {

    ## Write collapsed gene-level BED file to store flattened exons
    for(g in exons_g){
        for(r in exons_g[g]){
	    printf exons_g[g][r]"\n" >> "temp.bed";
	}
	system("bedtools sort -i temp.bed > sorted.bed");
	call="bedtools merge -s -o distinct -c 4,5,6 -i sorted.bed >> flat_exons.bed"
	system(call);
	close("temp.bed");
	close("sorted.bed");
	system("rm temp.bed");
	system("rm sorted.bed");
    }
    close("flat_exons.bed");
    
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
    flag=(getline line < "flat_exons.bed");
    idx=0;
    print(line)
    while(flag == 1){
	idx++;
	gene_bed[idx] = line;
	flag=(getline line < "flat_exons.bed");
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

    ## Calculate exon lengths and nucleotide content
    call="bedtools nuc -fi "genome" -bed flat_exons.bed";
    call=call" > nuc_exons.bed";
    system(call);

    ## Calculate Gene Level length and nucleotide counts
    call="bedtools groupby -i nuc_exons.bed ";
    call=call"-g 4 -o sum -c 9,10,11,12,13,14,15 > ";
    call=call"nuc_gene.tsv"
    system(call);

    ## Load gene level nucleotide counts
    line=""
    flag=(getline line < "nuc_gene.tsv");
    idx=0;
    while(flag == 1){
	idx++;
	nuc_gene[idx] = line;
	flag=(getline line < "nuc_gene.tsv");
    }
    close("nuc_gene.tsv")
    
    ## Calculate GC Percentage and write to file
    for(nuc_row in nuc_gene){
	m=split(nuc_gene[nuc_row],row,"\t");
	GC = (row[3] + row[4]) / row[8];
	new_row=row[1]"\t"row[8]"\t"GC;
	printf(new_row"\n") >> "Gene_Length_GC.tsv";
    }
    close("Gene_Length_GC.tsv");
    system("rm nuc_gene.tsv")
    system("rm nuc_exon.tsv")
}
