#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=2G
#$ -l h_rt=6:0:0
#$ -cwd
#$ -j y
#$ -N job-populations-together-all-snps-r80-default-1

module load use.dev
module add stacks/2.41

populations \
   -P /data/scratch/mpx469/STACKS/gstacks/gstacks-together-output/ \
   -O populations-together-all-snps-r80-default-1 \
   -M /data/scratch/mpx469/STACKS/gstacks/popmap-selection-together.txt \
   -r 0.8 \
   --vcf --fasta-loci --plink
   

