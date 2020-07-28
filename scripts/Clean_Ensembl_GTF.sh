#!/bin/bash
#
# File: Check_Ensembl_GTF.sh
# Author: Adam Faranda
# Purpose: Remove transcripts shorter than detectable from a reference GTF file
# Created: July 27, 2020
# 
# 
export igtf=Mus*.gtf  #Name of the input gtf file
export ogtf=out.gtf   #Name of the output gtf file
export nhln=5  #Number of header lines
export ml=1000  #Minimum number of nucleotides for a given transcript
awk -v NHLN=$nhln -v OFS="\t" -v ML=$ml -F"\t| |;"\
    '{gsub(/;\040+/,";")}\
	  (NR < NHLN)\
	       {next}
	  (NR == FNR && $3 == "exon")\
	      {if(tlen[$14] == "") tlen[$14] = ($5 - $4); else tlen[$10] += ($5 - $4)}\
	  (NR == FNR && $3 == "exon")\
	      {if(tlen[$14] < ML) a[$14]; 
	       else {\
	       	    delete a[$14];\
	       	    b[$10];\
	      	    {if(mn[$10] == ""){ mn[$10] = $4; mx[$10] = $5}};\
 	      	    {if($4 < mn[$10]) mn[$10] = $4};\
 	      	    {if($5 > mx[$10]) mx[$10] = $5}\
	       }}\
	  (NR == FNR)\
               {next}\
          (FNR < NHLN)\
 	       {next}\
          ($10 in b && $3 == "gene")\
	       {\
		$4=mn[$10];\
		$5=mx[$10];\
		gsub(";","; ");\
		for (i=1; i <= NF; i++){\
		    if(i < 9) printf $i"\t";\
		    else if(i > 8 && (i % 2) == 1) printf $i" ";\
		    else if(i > 8 && (i % 2) == 0) printf $i"; ";
		 }\
		 printf("\n");\
	       }\
	  ($3 != "gene" && !($14 in a))\
	       {gsub(";","; "); print $0}'\
	  $igtf $igtf | cut -d ";" -f1-5 > $ogtf


