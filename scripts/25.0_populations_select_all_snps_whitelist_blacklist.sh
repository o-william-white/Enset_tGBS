#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_populations_select_all_snps_whitelist_blacklist
#$ -t 70,80,90,100,110,120
#$ -tc 6

module load use.dev
module add stacks/2.41

# mkdir if not already present
mkdir -p populations/populations_select_${SGE_TASK_ID}_all_snps_whitelist_blacklist

populations \
   -P gstacks/gstacks_${SGE_TASK_ID} \
   -O populations/populations_select_${SGE_TASK_ID}_all_snps_whitelist_blacklist \
   -M gstacks/popmap_select_all.txt \
   -W populations/populations_select_${SGE_TASK_ID}_all_snps_whitelist/whitelist_overall \
   -B blacklist/select_overall_blasklist/blacklist_${SGE_TASK_ID}.txt \
   --vcf --fasta-loci --plink \

