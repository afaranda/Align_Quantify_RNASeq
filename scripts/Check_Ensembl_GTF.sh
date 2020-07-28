#!/bin/bash
#
# File: Check_Ensembl_GTF.sh
# Author: Adam Faranda
# Purpose: Verify that the input GTF file is properly formed
# Created: July 27, 2020
# 
# 
export igtf=out.gtf  #Name of the input gtf file
export nhln=1  #Number of header lines
export ml=1000  #Minimum number of nucleotides for a given transcript

# Verify that every subordinate feature has a gene id for a gene level feature
gawk -v NHLN=$nhln -v OFS="\t" -F"\t| |;"\
	  '{gsub(/;\040+/,";")}\
	  (NR > NHLN && $3 == "gene") {a[$10]; next} \
	  (NR > NHLN && $3 != "gene") {if(!($10 in a)) print $0; else next}'\
	  $igtf > feature_no_gene.txt

FNG=$(cat feature_no_gene.txt | wc -l | tr -d '[:space:]')
echo "Subordinate Features with no gene in raw GTF $FNG"

# Verify that columns 13 & 14 define a non-null
# transcript ID for all non-gene (subordinate) features
gawk -v NHLN=$nhln -v OFS="\t" -F"\t| |;"\
	  '{gsub(/;\040+/,";")}\
	  (NR > NHLN && $3 != "gene"){\
	       if($13 == "transcript_id" && $14 ~ /^\042ENSMUST/) next;\
	       else print $0}' $igtf > bad_txid.txt

BTX=$(cat bad_txid.txt | wc -l | tr -d '[:space:]')
echo "Subordinate Features with no transcript ID in raw GTF $BTX"

# Verify that all features for each gene are within that gene's bound
gawk -v NHLN=$nhln -v OFS="\t" -F"\t| |;"\
	  '{gsub(/;\040+/,";")}
	  (NR < NHLN){next}
	  ($3 == "gene")\
	  {
		if( !($10 in g))
		{
			g[$10];
			mng[$10]=$4;
			mxg[$10]=$5
		}
		else 
		{
			print $10" is a duplicated gene ID"
		}
	  }
	  ($3 == "transcript")\
	  {
		if( !($14 in t))
		{
			t[$14];
			mnt[$14]=$4;
			mxt[$14]=$5
		}
		else 
		{
			print $14" is a duplicated gene ID"
		}
	  } 

	  ($3 != "gene")\
	  {
		if($4 < mng[$10]) {print "bad feature " $10, $3, $4, $5}
		if($5 > mxg[$10]) {print "bad feature " $10, $3, $4, $5}
		if($3 == "transcript")
		{
		      if($4 < mnt[$14]) {print "bad feature " $10, $3, $4, $5}
		      if($5 > mxt[$14]) {print "bad feature " $10, $3, $4, $5}
			
		}
	  }' $igtf > "feature_out_of_bounds.txt"
FOB=$(cat feature_out_of_bounds.txt | wc -l | tr -d '[:space:]')
echo "Features outside Gene or Transcript bounds in raw GTF $FOB"


if [ $FNG == 0 ] && [ $BTX == 0 ] && [ $FOB == 0 ]; then
    echo "GTF Valid"
    rm bad_txid.txt feature_no_gene.txt 
else
    echo "Invalid GTF"
fi
