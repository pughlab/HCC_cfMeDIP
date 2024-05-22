# Description: Export wig files using MEDIPS.exportWIG() from the BAM files generated by step_1_fastq_to_bam.sh
# Navigate to the directory containing the wig files generated by bam_to_wig.R

cd /path/to/your/wig
mkdir ../matrices

echo " wig to matrices converting initiated......"

# Generate matrics
for i in *.wig

do

j=$i$".txt"

awk 'FNR > 2 {
    split($0, newArray, " ")
    if (length(newArray) > 1) print "0"  >> "'$j'"
    else if (length(newArray) == 1) print $0 >> "'$j'"
}' $i

echo "0" >> $j

echo $i$" done"

mv *.txt ../matrices

# Optional: Remove the original wig files to save space
# echo "Removing wig files..."
# rm *.wig

done
