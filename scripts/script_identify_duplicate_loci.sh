#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=6G
#$ -j y
#$ -l h_rt=1:0:0
#$ -N job_identifiy_duplicate_loci
#$ -t 1-12
#$ -tc 6

LENGTH=$(sed -n "${SGE_TASK_ID}p" input_args | cut -f 1 -d " " )
DATASET=$(sed -n "${SGE_TASK_ID}p" input_args | cut -f 2 -d " " )

module load R

Rscript identify_duplicate_loci.R /data/scratch/mpx469/tGBS_enset_project/populations/populations_${LENGTH}_${DATASET}_output/populations.snps.vcf duplicate_loci_${LENGTH}_${DATASET}_output/

echo done

