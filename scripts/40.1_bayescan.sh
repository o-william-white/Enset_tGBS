#!/bin/bash
#$ -pe smp 8
#$ -l h_vmem=1G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y
#$ -N job_bayescan
#$ -t 1-3

mkdir -p bayescan/bayescan_output_${SGE_TASK_ID}

/data/home/mpx469/software/bayescan/BayeScan2.1/source/bayescan_2.1 \
   bayescan/populations.snps.maf.bayescan.txt \
   -od bayescan/bayescan_output_${SGE_TASK_ID} \
   -threads 8

echo done

