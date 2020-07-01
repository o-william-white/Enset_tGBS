#!/bin/sh
#$ -cwd
#$ -j y
#$ -pe smp 24
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G
#$ -N job_gstacks_distant
#$ -t 70,80,90,100,110,120
#$ -tc 6

module load use.dev
module add stacks/2.41

gstacks \
   -I /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti/samtools_${SGE_TASK_ID}_output/ \
   -M popmap_distant.txt \
   -O gstacks_distant_${SGE_TASK_ID}_output \
   -t 12

