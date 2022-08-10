
# write popmap
awk -F "," ' { print $2 "\t" $1 } ' tGBS_metadata_phylogenetic_analysis.csv | tail -n +2 > gstacks/popmap.txt

