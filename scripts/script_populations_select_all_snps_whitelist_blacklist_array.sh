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
mkdir -p populations_select_${SGE_TASK_ID}_all_snps_whitelist_blacklist_output

populations \
   -P /data/scratch/mpx469/tGBS_enset_project/gstacks/gstacks_${SGE_TASK_ID}_output/ \
   -O populations_select_${SGE_TASK_ID}_all_snps_whitelist_blacklist_output \
   -M /data/scratch/mpx469/tGBS_enset_project/gstacks/popmap_select_all.txt \
   -W /data/scratch/mpx469/tGBS_enset_project/populations/populations_select_${SGE_TASK_ID}_all_snps_output/whitelist_overall \
   -B /data/scratch/mpx469/tGBS_enset_project/blacklists/blacklists_overall/blacklist_select_${SGE_TASK_ID}_all_snps.txt \
   --vcf --fasta-loci --plink \

