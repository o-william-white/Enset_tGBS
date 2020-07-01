#!/bin/bash
#$ -pe smp 8
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_bayescan_example

mkdir -p bayescan_output

/data/home/mpx469/software/bayescan/BayeScan2.1/source/bayescan_2.1 test_SNPs.txt -od bayescan_output -threads 8

echo done

