#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_populations_distant_single_snp_blacklist
#$ -t 70,80,90,100,110,120
#$ -tc 6

module load use.dev
module add stacks/2.41

populations \
   -P /data/scratch/mpx469/tGBS_enset_project/gstacks/gstacks_distant_${SGE_TASK_ID}_output/ \
   -O populations_distant_${SGE_TASK_ID}_single_snp_blacklist_output \
   -M /data/scratch/mpx469/tGBS_enset_project/gstacks/popmap_distant.txt \
   -R 0.8 \
   --min-mac 3 \
   -B /data/scratch/mpx469/tGBS_enset_project/blacklists/blacklist_overall_distant_${SGE_TASK_ID}_single_snp \
   --write-single-snp \
   --phylip --phylip-var --phylip-var-all --vcf --fasta-loci --plink --genepop
   