#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_estimate_mlg

ml R/3.6.1

# estimate mlg
Rscript additional_scripts/estimate_mlg.R

# check that mlgs are monophyletic
Rscript additional_scripts/check_mlg_monophyly.R

echo Complete!

