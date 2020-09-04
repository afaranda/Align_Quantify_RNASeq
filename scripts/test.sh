# echo "'"
# while IFS= read -r line
# do
#     echo "(\$0 ~ "${line}") {if(\$0 ~ \"transcript_id\") print \$0;"\
#     "else print \$0\" transcript_id \"\";\"}"\\
# done < $1
gawk -v FS="\t"\
     'BEGIN {rc=0; fc=0}
     (rc > FNR)\
     {
        fc ++
        if(fc == 1)
        {
           unl[$1]
        }
        else
        {
           rib[$1]
        }
     }
    {rc = FNR - 1}
    (fc == 0){aln[$1]}
    (fc == 1){unl[$1]}
    (fc == 2){rib[$1]}

    END\
    {
        alnCount = 0
        unlCount = 0
        ribCount = 0
        nulCount = 0

        for(read in rib)
        {
            ribCount ++
            if(read in aln)
            {
               alnCount ++
            }
            else if(read in unl)
            {
               unlCount ++
            }
            else
            {
               nulCount ++
            }
        }
        print ribCount "\t" alnCount "\t" unlCount "\t" nulCount "\t"
    }' $1 $2 $3 
