#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=04:0:0
#$ -l h_vmem=1G
#$ -N job-read-number-count

# count raw read number from data2bio
cat /data/scratch/mpx469/sample-list.txt | while read i; do
   zcat /data/scratch/mpx469/Data2Bio_final/raw/$i'.digested.fq.gz' | echo $((`wc -l`/4)) >> count-data2bio-raw.txt
done

# count trimmed read number from data2bio
cat /data/scratch/mpx469/sample-list.txt | while read i; do
   zcat /data/scratch/mpx469/Data2Bio_final/trimmed/$i'.digested.trimmed.fq.gz' | echo $((`wc -l`/4)) >> count-data2bio-trimmed.txt
done

# count read number after running trimmomatic
cat /data/scratch/mpx469/sample-list.txt | while read i; do
   zcat /data/scratch/mpx469/trimmomatic/trimmomatic-output/$i'.digested.trimmomatic.fq.gz' | echo $((`wc -l`/4)) >> count-trimmomatic.txt
done

