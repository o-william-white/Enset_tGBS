#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_phylosignal

ml R/3.6.1

mkdir -p phylosignal

# phylosignal
Rscript additional_scripts/phylosignal.R

echo Complete!

