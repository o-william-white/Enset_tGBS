#!/bin/sh
#$ -cwd
#$ -j y
#$ -l h_rt=2:0:0
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -t 1-283
#$ -tc 283
#$ -N job-bwa-map

module add bwa

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" /data/scratch/mpx469/sample-list.txt)

cp /data/scratch/mpx469/trimmomatic/trimmomatic-output/${INPUT_FILE}.digested.trimmomatic.fq.gz ./trimmomatic-fq/

gunzip ./trimmomatic-fq/${INPUT_FILE}.digested.trimmomatic.fq.gz

bwa mem GCA_000331365.3_Ensete_JungleSeeds_v3.0_genomic.fna ./trimmomatic-fq/${INPUT_FILE}.digested.trimmomatic.fq > bwa-map-output/${INPUT_FILE}.mapped.sam

