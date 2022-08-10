#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_populations_select_all_snps_whitelist
#$ -t 70,80,90,100,110,120
#$ -tc 6

module load use.dev
module add stacks/2.41

# mkdir if not already present
mkdir -p populations/populations_select_${SGE_TASK_ID}_all_snps_whitelist

# create overall whitelist
cat populations/populations_select_${SGE_TASK_ID}_all_snps/all/whitelist \
    populations/populations_select_${SGE_TASK_ID}_all_snps/dom/whitelist \
    populations/populations_select_${SGE_TASK_ID}_all_snps/wil/whitelist | grep -e "Locus ID" -v | sort -n | uniq > populations/populations_select_${SGE_TASK_ID}_all_snps_whitelist/whitelist_overall

# populations using whitelist
populations \
   -P gstacks/gstacks_${SGE_TASK_ID} \
   -O populations/populations_select_${SGE_TASK_ID}_all_snps_whitelist \
   -M gstacks/popmap_select_all.txt \
   -W populations/populations_select_${SGE_TASK_ID}_all_snps_whitelist/whitelist_overall \
   --vcf --fasta-loci 

echo Complete!

