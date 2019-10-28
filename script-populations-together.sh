#!/bin/bash
#$ -pe smp 12
#$ -l h_vmem=4G
#$ -l h_rt=2:0:0
#$ -cwd
#$ -j y
#$ -N job-populations-together

module load use.dev
module add stacks/2.41

populations \
   -P /data/scratch/mpx469/stacks/ref-map/gstacks/gstacks-together-output/ \
   -O populations-together-output \
   -M /data/scratch/mpx469/stacks/ref-map/gstacks/popmap-selection-together.txt \
   -r 0.4 \
   --vcf --plink --structure \
   -t 12

