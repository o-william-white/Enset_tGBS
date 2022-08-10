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
mkdir -p blacklist/select_high_depth/
mkdir -p blacklist/select_high_depth/high_depth_${SGE_TASK_ID}

Rscript additional_scripts/high_depth_loci.R populations/populations_select_${SGE_TASK_ID}_all_snps_whitelist/populations.snps.vcf blacklist/select_high_depth/high_depth_${SGE_TASK_ID}

echo done

