#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_run_top_go

ml R/3.6.1

Rscript additional_scripts/run_topGO.R

echo Complete!


