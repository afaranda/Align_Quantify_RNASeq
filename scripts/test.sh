echo "'"
while IFS= read -r line
do
    echo "(\$0 ~ "${line}") {if(\$0 ~ \"transcript_id\") print \$0;"\
    "else print \$0\" transcript_id \"\";\"}"\\
done < $1

