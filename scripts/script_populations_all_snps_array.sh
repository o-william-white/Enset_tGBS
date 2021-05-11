#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_populations_all_snps
#$ -t 70,80,90,100,110,120
#$ -tc 6

module load use.dev
module add stacks/2.41

# mkdir if not already present
mkdir -p populations_${SGE_TASK_ID}_all_snps_output

populations \
   -P /data/scratch/mpx469/tGBS_enset_project/gstacks/gstacks_${SGE_TASK_ID}_output/ \
   -O populations_${SGE_TASK_ID}_all_snps_output \
   -M /data/scratch/mpx469/tGBS_enset_project/gstacks/popmap.txt \
   -R 0.8 \
   --min-mac 3 \
   --vcf --fasta-loci
   

