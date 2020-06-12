#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=2:0:0
#$ -l h_vmem=1G
#$ -j y
#$ -N job_00_count_raw_reads

# count read number after running raw
cat /data/scratch/mpx469/tGBS_enset_project/sample_list.txt | while read i; do
   zcat /data/scratch/mpx469/tGBS_enset_project/Data2Bio_final/raw/${i}.digested.fq.gz | echo $((`wc -l`/4)) >> count_00_raw.txt
done

echo done

