#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_blastn_musa

ml blast+
ml R/3.6.1

# downloaded from banana genome hub for Musa_acuminata_ssp_malaccensis_2.0
# https://banana-genome-hub.southgreen.fr/node/50/7720981
cp additional_data/musa_acuminata_v2_cds.fna bedtools/

# makeblastdb
makeblastdb -in bedtools/outliers_10kb.fasta   -dbtype nucl
makeblastdb -in bedtools/background_10kb.fasta -dbtype nucl

# blast
blastn -query bedtools/musa_acuminata_v2_cds.fna -db bedtools/outliers_10kb.fasta   -outfmt 6 -evalue 1e-30 > bedtools/blast_out_outliers
blastn -query bedtools/musa_acuminata_v2_cds.fna -db bedtools/background_10kb.fasta -outfmt 6 -evalue 1e-30 > bedtools/blast_out_background

# get top hits 
Rscript additional_scripts/top_hit.R bedtools/blast_out_outliers
Rscript additional_scripts/top_hit.R bedtools/blast_out_background

echo Complete!

