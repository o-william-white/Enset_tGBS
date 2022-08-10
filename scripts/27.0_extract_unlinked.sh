#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_extract_unlinked_select
#$ -t 1-6

module load anaconda3
source activate stacks_workflow
ml R/3.6.1

SNP_TYPE=$(echo -e "diverged\nduplicated\nhighcov\nlowconf\nmas\nsingleton" | sed -n "${SGE_TASK_ID}p")

mkdir -p extract_unlinked_select

# extract unlinked
# minimum difference between 0.0 and 1.0 to keep a second snp = 0.5 (recommended)
# maximum distance in base pairs to consider linked SNPs = 10kb 
python3 \
  additional_scripts/stacks_workflow/00-scripts/11_extract_unlinked_snps_genome.py \
  snp_duplication_select/populations.snps.${SNP_TYPE}.vcf \
  0.5 \
  10000 \
  extract_unlinked_select/populations.snps.${SNP_TYPE}.unlinked.vcf

# count the number of singleton snps
N_SING=$(grep -e "^#" -v -c snp_duplication_select/populations.snps.${SNP_TYPE}.vcf)

# count the number of snps retained
N_SNPS=$(grep -e "^#" -v -c extract_unlinked_select/populations.snps.${SNP_TYPE}.unlinked.vcf)

# print the number of snps
echo ${SNP_TYPE}
echo ${N_SING} snps
echo ${N_SNPS} unlinked snps retained

echo done

