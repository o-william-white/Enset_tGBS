#!/bin/sh
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -t 1-283
#$ -tc 36
#$ -N job_cutadapt_array

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" /data/scratch/mpx469/tGBS_enset_project/sample_list.txt)

source /data/home/mpx469/software/python-virtualenv/cutadapt/bin/activate

cutadapt \
   -g ^CATG \
   -g ^GATC \
   --action=none \
   --untrimmed-output cutadapt_output/${INPUT_FILE}.cutadapt.without.cutsites.fq.gz \
   -o cutadapt_output/${INPUT_FILE}.cutadapt.with.cutsites.fq.gz \
   /data/scratch/mpx469/tGBS_enset_project/trimmomatic/trimmomatic_output/${INPUT_FILE}.trimmomatic.fq.gz

echo done

