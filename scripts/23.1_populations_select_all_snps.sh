#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_populations_select_all_snps
#$ -t 1-18

module load use.dev
module add stacks/2.41

LENGTH=$(sed -n "${SGE_TASK_ID}p" input_args_select | cut -f 1)
POPMAP=$(sed -n "${SGE_TASK_ID}p" input_args_select | cut -f 2)

# mkdir if not already present
mkdir -p populations/populations_select_${LENGTH}_all_snps/
mkdir -p populations/populations_select_${LENGTH}_all_snps/${POPMAP}

populations \
   -P gstacks/gstacks_${LENGTH} \
   -O populations/populations_select_${LENGTH}_all_snps/${POPMAP} \
   -M gstacks/popmap_select_${POPMAP}.txt \
   -R 0.8 \
   --vcf --fasta-loci

# write whitelist
grep -e "^# pop" -v populations/populations_select_${LENGTH}_all_snps/${POPMAP}/populations.sumstats.tsv | cut -f 1,4 | sed 's/# Locus ID/Locus ID/g' | sort -n | uniq > populations/populations_select_${LENGTH}_all_snps/${POPMAP}/whitelist

echo Complete!

