#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=01:0:0
#$ -l h_vmem=1G
#$ -j y
#$ -N job_05_count_samtools
#$ -t 70,80,90,100,110,120

module load samtools

# count read number after running bamtools
cat /data/scratch/mpx469/tGBS_enset_project/sample_list.txt | while read i; do

   # count unmapped reads
   samtools view -f 0x4 /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti/samtools_${SGE_TASK_ID}_output/${i}.unique.sorted.bam  | wc -l >> count_05_samtools_${SGE_TASK_ID}_unmapped.txt

   # count mapped reads
   samtools view -F 0x4 /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti/samtools_${SGE_TASK_ID}_output/${i}.unique.sorted.bam  | wc -l >> count_05_samtools_${SGE_TASK_ID}_mapped.txt

done

echo done

