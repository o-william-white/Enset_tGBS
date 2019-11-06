#!/bin/bash
#$ -pe smp 12
#$ -l h_vmem=24G
#$ -l h_rt=10:0:0
#$ -l node_type=sm
#$ -l highmem
#$ -cwd
#$ -j y
#$ -N job-populations-separate

module load use.dev
module add stacks/2.41

populations \
   -P /data/scratch/mpx469/stacks/ref-map/gstacks/gstacks-separate-output/ \
   -O populations-separate-output \
   -M /data/scratch/mpx469/stacks/ref-map/gstacks/popmap-selection-separate.txt \
   -R 0.8 \
   --write-single-snp \
   --phylip --phylip-var --phylip-var-all --vcf \
   -t 12

