#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_format_top_go_input

mkdir -p top_go

# downloaded from banana genome hub for Musa_acuminata_ssp_malaccensis_2.0
# https://banana-genome-hub.southgreen.fr/node/50/7720981 
#cp additional_data/musa_acuminata_v2.gff3 top_go

# format GO annotations 
# http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html
#grep "GO:" top_go/musa_acuminata_v2.gff3 | cut -f 9 | while read LINE; do
#  ID=$(echo $LINE | cut -f 1 -d ";" | sed 's/ID=//g')
#  GO=$(echo $LINE | cut -f 6 -d ";" | sed 's/ISS_[0-9],GO:/\tGO:/g' | cut -f 2 | sed 's/,/, /g' )
#  echo -e "$ID\t$GO" 
#done | sort -k 1,1 -V > top_go/annotations_ID2GO.txt

# genes of interest
# tail -n +2 bedtools/blast_out_outliers_top_hits | cut -f 1 > top_go/genes_of_interest.txt

# background genes
tail -n +2 bedtools/blast_out_background_top_hits | cut -f 1 | grep -f - top_go/annotations_ID2GO.txt > top_go/annotations_background.txt

echo Complete!

