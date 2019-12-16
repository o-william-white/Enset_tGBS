#!/bin/sh
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -t 1-283
#$ -c 36
#$ -N job-trimmomatic-array

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" /data/scratch/mpx469/sample-list.txt)

module load trimmomatic/0.36

java -jar /share/apps/centos7/trimmomatic/0.36/trimmomatic-0.36.jar  \
   SE \
   /data/scratch/mpx469/Data2Bio_final/raw/${INPUT_FILE}.digested.fq.gz \
   /data/scratch/mpx469/trimmomatic/trimmomatic-output/${INPUT_FILE}.trimmomatic.fq.gz \
   LEADING:15 \
   TRAILING:15 \
   SLIDINGWINDOW:4:15 \
   MINLEN:36


