#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_identify_sites_with_high_depth_array
#$ -t 1-12

LENGTH=$(sed -n "${SGE_TASK_ID}p" input_args | cut -f 1 -d " ")
DATASET=$(sed -n "${SGE_TASK_ID}p" input_args | cut -f 2 -d " ")

module load R

Rscript identify_sites_with_high_depth.R ${LENGTH} ${DATASET}

echo done

