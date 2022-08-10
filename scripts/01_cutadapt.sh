#!/bin/sh
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -t 1-283
#$ -tc 36
#$ -N job_cutadapt

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" sample_list.txt)

source /data/home/mpx469/software/python-virtualenv/cutadapt/bin/activate

# mkdir if not already present
mkdir -p cutadapt

cutadapt \
   -g ^CATG \
   -g ^GATC \
   --action=none \
   --untrimmed-output cutadapt/${INPUT_FILE}.without.cutsites.fq.gz \
   -o cutadapt/${INPUT_FILE}.with.cutsites.fq.gz \
   trimmomatic/${INPUT_FILE}.fq.gz

echo done

