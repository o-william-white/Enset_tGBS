#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=2:0:0
#$ -l h_vmem=1G
#$ -j y
#$ -N job_03_count_process_radtags_reads
#$ -t 70,80,90,100,110,120

# count read number after running process_radtags
cat /data/scratch/mpx469/tGBS_enset_project/sample_list.txt | while read i; do
   zcat /data/scratch/mpx469/tGBS_enset_project/process_radtags/process_radtags_${SGE_TASK_ID}_output/${i}.cutadapt.with.cutsites.fq.gz | echo $((`wc -l`/4)) >> count_03_process_radtags_${SGE_TASK_ID}.txt
done

echo done

