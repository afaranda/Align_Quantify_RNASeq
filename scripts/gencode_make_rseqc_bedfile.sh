#!/bin/bash

GTFPATH=/work/abf/MouseGencodeM25/gencode.vM25.basic.annotation.gtf
BEDPATH=/work/abf/MouseGencodeM25/rseqc_gene_models.bed

# Convert selected genes from GTF File to bed to get gene models for RSeQC
# Add genes as rows: ($0 ~ XXXX) where XXXXX is an ensembl gene ID
gawk '($0 ~ "ENSMUSG00000027168")\
    {if($0 ~ "transcript_id") print $0;\
    else print $0" transcript_id \"\";"}' ${GTFPATH}\
    | gtf2bed - > converted.bed

gawk -v OFS="\t" -F"\t| |;"\
     'BEGIN{first_tx=0}
      {gsub("; ",";")}\
      (NR == FNR && $8 != "exon")\
      	  {next}\
      (NR == FNR && $8 == "exon")\
          {ex_start[$13][$25]=$2; ex_end[$13][$25]=$3;next}\
      (NR != FNR && $8 == "transcript")\
	  {  	
		exc=0
		m=asort(ex_start[$13],exs,"@val_num_asc")
		n=asort(ex_end[$13],exe,"@val_num_asc")
		gsub("\042","",$13)
		if(first_tx == 0){
	  	       printf($1"\t"$2"\t"$3"\t"$13"\t""0""\t"$6"\t"$2"\t"$3"\t""0""\t"m"\t")
		       for(i=1; i <=m; i++){
		       	     printf((exe[i] - exs[i])",")
		       }
		       printf("\t")
		       for(i=1; i <=m; i++){
		       	     printf((exs[i] - $2)",")
		       }
		       first_tx=1
		}
		else{
	  	       printf("\n"$1"\t"$2"\t"$3"\t"$13"\t""0""\t"$6"\t"$2"\t"$3"\t""0""\t"m"\t")
		       for(i=1; i <=m; i++){
		       	     printf((exe[i] - exs[i])",")
		       }
		       printf("\t")
		       for(i=1; i <=m; i++){
		       	     printf((exs[i] - $2)",")
		       }
		}
		$13 = "\042" $13 "\042"
	  }' converted.bed converted.bed > $BEDPATH
rm converted.bed
