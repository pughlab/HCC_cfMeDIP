# Description: This script extract the ISD (insert size distribution) from fastp report which can be used to generate Fig.S2B & 2C
# Author：Kui Chen → kui.chen@uhn.ca  ###########################
# 1st Created: 2021-10-01  #####################################
# Last Modified: 2023-12-10  ###################################


#!/bin/bash

files=(*.html)
tmp_files=()

for file in "${files[@]}"
do
    prefix=$(echo $file | sed 's/.html//')
    output_file="${prefix}_output.csv"
    tmp_files+=($output_file)

    # Prepare the header with custom column names
    header="${prefix}_size(bp),${prefix}_%"

    # Extract and process data, prepend header
    echo $header > $output_file
    awk '/plot_insert_size/{count++} count==1 && $0 !~ /plot_insert_size/{print} count>1{exit}' $file | grep 'var data=' |
    sed -n 's/.*x:\[\(.*\)\],y:\[\(.*\)\].*/\1 \2/p' |
    awk -v prefix=$prefix '{
        split($1, x, ",");
        split($2, y, ",");
        for (i=1; i<=length(x); i++) print x[i] ", " y[i]
    }' >> $output_file  # Note the append redirection '>>' to keep the header

done

# Combine all output files by columns
paste -d ',' ${tmp_files[@]} > ISD.csv

# Clean up temporary files
# Uncomment the next line to delete individual output files after combining
rm ${tmp_files[@]}

echo "Data extraction and merging complete. Combined output file: combined_output.csv"
