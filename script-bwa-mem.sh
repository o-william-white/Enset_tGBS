#!/bin/sh
#$ -cwd
#$ -j y
#$ -l h_rt=2:0:0
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -t 1-283
#$ -tc 36
#$ -N job-bwa-mem

module add bwa

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" /data/scratch/mpx469/sample-list.txt)

cp /data/scratch/mpx469/cutadapt/cutadapt-output/${INPUT_FILE}.cutadapt.with.cutsites.fq.gz ./cutadapt-fq/

gunzip ./cutadapt-fq/${INPUT_FILE}.cutadapt.with.cutsites.fq.gz

bwa mem GCA_000331365.3_Ensete_JungleSeeds_v3.0_genomic.fna ./cutadapt-fq/${INPUT_FILE}.cutadapt.with.cutsites.fq > bwa-mem-output/${INPUT_FILE}.mapped.sam

