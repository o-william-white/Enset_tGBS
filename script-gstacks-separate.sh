#!/bin/sh
#$ -cwd
#$ -j y
#$ -pe smp 12
#$ -l h_rt=12:0:0
#$ -l h_vmem=2G
#$ -N job-gstacks-separate

module load use.dev
module add stacks/2.41

gstacks \
   -I /data/scratch/mpx469/stacks/ref-map/samtools/samtools-output/ \
   -M popmap-selection-separate.txt \
   -O gstacks-separate-output \
   -t 12

