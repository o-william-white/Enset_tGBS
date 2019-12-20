#!/bin/sh
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G
#$ -N job-summary-flagstat-counts

cat /data/scratch/mpx469/sample-list.txt | while read i; do
   head -n 1 samtools-output/${i}.flagstat.txt | cut -f 1 -d " ";
done > flagstat-bwa-total.txt

cat /data/scratch/mpx469/sample-list.txt | while read i; do
   awk 'NR==5' samtools-output/${i}.flagstat.txt | cut -f 1 -d " ";
done > flagstat-bwa-mapped.txt

cat /data/scratch/mpx469/sample-list.txt | while read i; do
   head -n 1 samtools-output/${i}.flagstat.filtered.txt | cut -f 1 -d " ";
done > flagstat-samtools-total.txt

cat /data/scratch/mpx469/sample-list.txt | while read i; do  
   awk 'NR==5' samtools-output/${i}.flagstat.filtered.txt | cut -f 1 -d " "; 
done > flagstat-samtools-mapped.txt

