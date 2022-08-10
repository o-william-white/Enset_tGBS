#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=6G
#$ -j y
#$ -l h_rt=1:0:0
#$ -N job_high_depth_loci
#$ -t 70,80,90,100,110,120

ml R/3.6.1

# mkdir if not already present
mkdir -p blacklist/high_depth/
mkdir -p blacklist/high_depth/high_depth_${SGE_TASK_ID}

Rscript additional_scripts/high_depth_loci.R /data/scratch/mpx469/tGBS_enset_project/populations/populations_${SGE_TASK_ID}_all_snps/populations.snps.vcf blacklist/high_depth/high_depth_${SGE_TASK_ID}

echo done

