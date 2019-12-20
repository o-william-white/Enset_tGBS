#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=2G
#$ -l h_rt=6:0:0
#$ -cwd
#$ -j y
#$ -N job-populations-separate-all-snps-R80-maf-het-2-blacklist

module load use.dev
module add stacks/2.41

populations \
   -P /data/scratch/mpx469/STACKS/gstacks/gstacks-separate-output/ \
   -O populations-separate-all-snps-R80-maf-het-2-blacklist \
   -M /data/scratch/mpx469/STACKS/gstacks/popmap-selection-separate.txt \
   -R 0.8 \
   -B /data/scratch/mpx469/STACKS/blacklist-duplicates/blacklist-separate-all-snps-R80-maf-het/blacklist.txt \
   --min-maf 0.05 \
   --max-obs-het 0.8 \
   --phylip --phylip-var --phylip-var-all --vcf --fasta-loci --plink
   

