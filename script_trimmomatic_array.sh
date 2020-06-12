#!/bin/sh
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -t 1-283
#$ -tc 36
#$ -N job_trimmomatic_array

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" /data/scratch/mpx469/tGBS_enset_project/sample_list.txt)

module load trimmomatic/0.36

java -jar /share/apps/centos7/trimmomatic/0.36/trimmomatic-0.36.jar  \
   SE \
   /data/scratch/mpx469/tGBS_enset_project/Data2Bio_final/raw/${INPUT_FILE}.digested.fq.gz \
   /data/scratch/mpx469/tGBS_enset_project/trimmomatic/trimmomatic_output/${INPUT_FILE}.trimmomatic.fq.gz \
   LEADING:15 \
   TRAILING:15 \
   SLIDINGWINDOW:4:15 \
   MINLEN:36

echo done

