#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=2:0:0
#$ -l h_vmem=1G
#$ -j y
#$ -N job_01_count_trimmomatic_reads

# count read number after running trimmomatic
cat /data/scratch/mpx469/tGBS_enset_project/sample_list.txt | while read i; do
   zcat /data/scratch/mpx469/tGBS_enset_project/trimmomatic/trimmomatic_output/${i}.trimmomatic.fq.gz | echo $((`wc -l`/4)) >> count_01_trimmomatic.txt
done

echo done

