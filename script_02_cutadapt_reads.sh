#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=2:0:0
#$ -l h_vmem=1G
#$ -j y
#$ -N job_02_count_cutadapt_reads

# count read number after running cutadapt
cat /data/scratch/mpx469/tGBS_enset_project/sample_list.txt | while read i; do
   zcat /data/scratch/mpx469/tGBS_enset_project/cutadapt/cutadapt_output/${i}.cutadapt.with.cutsites.fq.gz | echo $((`wc -l`/4)) >> count_02_cutadapt.txt
done

echo done

