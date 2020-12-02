#!/bin/bash
#SBATCH --job-name=make_rseqc_bed
#SBATCH --mem=32000
#SBATCH --ntasks=4
#
# File: targets_from_gtf.sh
# Purpose: Extract lines from a GTF file by matching an attribute
# Created: August 25, 2020
# Author: Adam Faranda
#
# Set Environment Variables
#  TARGETS= << Path to file listing Target ID's or Symbols in single column >>
#  TGTLABL= << Label appended to target ID, eg gene_id, gene_name etc. >>
#  GTFPATH= << Lines extracted from the file at this path >>
#  RESULTS= << Output stored in the file at this path >>
#

TARGETS=/work/abf/MouseEnsembl101/Crystallin_symbols.txt
TGTLABL=gene_name
GTFPATH=/work/abf/MouseEnsembl101/Mus_musculus.GRCm38.101.gtf
RESULTS=/work/abf/MouseEnsembl101/crystallins.gtf

if [ ! -z $TGTLABL ] && [ ! -z $TARGETS ]
then
    gawk -v lbl=${TGTLABL}\
	 -v FS="\t| |;"\
	 '\
          # Store values from first file in an array of label target pairs
          # with the target surrounded by quotes using \042 
	  (NR == FNR)\
          {
              tgt[lbl" \042"$1"\042"];
              next
          }
          
          # For lines from second (gtf) file iterate over array of targets
          # and print the line if a match exists
          (NR != FNR)\
          {
	     for(i in tgt){
                if($0 ~ i){
	  	    print $0
	  	}
             }
          }
          ' $TARGETS $GTFPATH > ${RESULTS}
fi
