#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=6G
#$ -j y
#$ -l h_rt=1:0:0
#$ -N job_duplicate_loci_select
#$ -t 70,80,90,100,110,120

ml R/3.6.1

mkdir -p blacklist/select_duplicate_loci
mkdir -p blacklist/select_duplicate_loci/duplicate_loci_${SGE_TASK_ID}

Rscript additional_scripts/duplicate_loci.R populations/populations_select_${SGE_TASK_ID}_all_snps_whitelist/populations.snps.vcf blacklist/select_duplicate_loci/duplicate_loci_${SGE_TASK_ID}/blacklist.txt

echo done

