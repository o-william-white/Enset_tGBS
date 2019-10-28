#!/bin/sh
#$ -cwd
#$ -j y
#$ -l h_rt=1:0:0
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -N job-samtools-flagstat
#$ -t 1-283
#$ -tc 283

module add samtools

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" /data/scratch/mpx469/sample-list.txt)

samtools flagstat /data/scratch/mpx469/stacks/ref-map/bwa/bwa-map-output/${INPUT_FILE}.mapped.sam > samtools-output/${INPUT_FILE}.flagstat.txt

samtools flagstat /data/scratch/mpx469/stacks/ref-map/samtools/samtools-output/${INPUT_FILE}.mapped.unique.sam > samtools-output/${INPUT_FILE}.flagstat.filtered.txt

