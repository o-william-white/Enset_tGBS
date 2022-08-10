#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_extrapolate_mlg

ml R/3.6.1

# estimate accumulation curve with rarefaction
Rscript additional_scripts/estimate_mlg_accumulation_curve.R

# extrapolate accumulation curve
# Rscript additional_scripts/estimate_mlg_accumulation_curve_extrapolate.R

echo Complete!

