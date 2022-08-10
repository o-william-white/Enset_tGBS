#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_dsuite_by_sample_results

ml R/3.6.1

Rscript additional_scripts/dsuite_by_sample_results.R

echo Complete!

