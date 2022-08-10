#!/bin/sh
#$ -cwd
#$ -j y
#$ -pe smp 12
#$ -l h_rt=1:0:0
#$ -l h_vmem=2G
#$ -N job_gstacks
#$ -t 70,80,90,100,110,120
#$ -tc 6

module load stacks/2.41

# mkdir if not already present
mkdir -p gstacks
mkdir -p gstacks/gstacks_${SGE_TASK_ID}

gstacks \
   -I map_reads_bedadeti/samtools_${SGE_TASK_ID}/ \
   -M gstacks/popmap.txt \
   -O gstacks/gstacks_${SGE_TASK_ID} \
   -t 12

