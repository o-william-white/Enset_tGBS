#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -j y
#$ -l h_rt=1:0:0
#$ -N job_write_overall_blacklist
#$ -t 70,80,90,100,110,120

ml R/3.6.1

mkdir -p blacklist/select_overall_blasklist

Rscript additional_scripts/write_overall_blacklist_select.R ${SGE_TASK_ID}

echo done

