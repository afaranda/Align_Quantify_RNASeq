#!/bin/bash

PREFIX=EnsMm_104

hisat2-inspect -s $PREFIX \
    | gawk -vOFS="\t" '($1 ~ /^Seq.*/) {print $2, $6}' > ${PREFIX}_Chrom_Sizes.txt

	 
