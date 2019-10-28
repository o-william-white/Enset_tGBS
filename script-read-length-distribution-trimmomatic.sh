#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=2:0:0
#$ -cwd
#$ -j y
#$ -N job-read-length-distribution-trimmomatic

cat /data/scratch/mpx469/sample-list.txt | while read i; do
   zcat /data/scratch/mpx469/trimmomatic/trimmomatic-output/$i'.digested.trimmomatic.fq.gz' | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > output-read-length-distribution-trimmomatic/read_length.$i'.digested.trimmomatic.txt'
done

