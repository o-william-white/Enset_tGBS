#!/bin/sh
#$ -cwd
#$ -j y
#$ -pe smp 12
#$ -l h_rt=12:0:0
#$ -l h_vmem=2G
#$ -N job-gstacks-together

module load use.dev
module add stacks/2.41

gstacks \
   -I /data/scratch/mpx469/stacks/ref-map/samtools/samtools-output/ \
   -M popmap-selection-together.txt \
   -O gstacks-together-output \
   -t 12

