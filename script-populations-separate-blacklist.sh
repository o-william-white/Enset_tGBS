#!/bin/bash
#$ -pe smp 12
#$ -l h_vmem=4G
#$ -l h_rt=10:0:0
#$ -cwd
#$ -j y
#$ -N job-populations-separate-blacklist

module load use.dev
module add stacks/2.41

populations \
   -P /data/scratch/mpx469/stacks/ref-map/gstacks/gstacks-separate-output/ \
   -O populations-separate-blacklist-output \
   -M /data/scratch/mpx469/stacks/ref-map/gstacks/popmap-selection-separate.txt \
   -R 0.8 \
   -B /data/scratch/mpx469/stacks/ref-map/filter-duplicates-separate/blacklist.txt \
   --write-single-snp \
   --phylip --phylip-var --phylip-var-all --vcf \
   -t 12

