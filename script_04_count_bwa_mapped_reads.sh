#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=01:0:0
#$ -l h_vmem=1G
#$ -j y
#$ -N job_04_count_bwa
#$ -t 70,80,90,100,110,120

module load samtools

# count read number after running bwa
cat /data/scratch/mpx469/tGBS_enset_project/sample_list.txt | while read i; do

   # count unmapped reads
   samtools view -f 0x4 /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti/bwa_mem_${SGE_TASK_ID}_output/${i}.bam  | wc -l >> count_04_bwa_mem_${SGE_TASK_ID}_unmapped.txt

   # count mapped reads
   samtools view -F 0x4 /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti/bwa_mem_${SGE_TASK_ID}_output/${i}.bam  | wc -l >> count_04_bwa_mem_${SGE_TASK_ID}_mapped.txt

   # count XA (alternative hits) and SA (chimeric) mapped reads
   samtools view /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti/bwa_mem_${SGE_TASK_ID}_output/${i}.bam | grep -e 'XA:Z:' -e 'SA:Z:' -c >> count_04_bwa_mem_${SGE_TASK_ID}_XA_SA.txt

done

echo done


