#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=2G
#$ -l h_rt=6:0:0
#$ -cwd
#$ -j y
#$ -N job-populations-separate-single-snp-R80-default-1

module load use.dev
module add stacks/2.41

populations \
   -P /data/scratch/mpx469/STACKS/gstacks/gstacks-separate-output/ \
   -O populations-separate-single-snp-R80-default-1 \
   -M /data/scratch/mpx469/STACKS/gstacks/popmap-selection-separate.txt \
   -R 0.8 \
   --write-single-snp \
   --phylip --phylip-var --phylip-var-all --vcf --fasta-loci --plink
   

