#!/bin/bash
#SBATCH --job-name=make_rseqc_bed
#SBATCH --mem=32000
#SBATCH --ntasks=4

TARGETS=/work/abf/MouseEnsembl100/Crystallin_symbols.txt
TGTLABL=gene_name
GTFPATH=/work/abf/MouseEnsembl100/cryaa.gtf

if [ ! -z $TGTLABL ] && [ ! -z $TARGETS ]
then
    echo cool
    gawk -v lbl=${TGTLABL}\
	 -v FS="\t| |;"\
	 '(NR == FNR) {tgt[$1]; next}
          (NR != FNR)\
              {
                 for(i=0; i<=NF; i++){                     
                     if($i == lbl){
                         #print($i " " $(i+1))
                         gsub("\042","",$(i+1))
                         if($(i+1) in tgt){
                             $(i+1)="\042"$(i+1)"\042"
                             for(j=1; j<=NF;j++){
                                if(j < 9) {
                                    printf($j"\t")
                                }
                                else if( j % 2 > 0){
                                    printf($j" ")
                                }
                                else if( j != NF){
                                    printf($j"; ")
                                }
                                else {
                                    print($j)
                                }
                             }
                         
                         }
                     }
                 }
              }' $TARGETS $GTFPATH
fi
