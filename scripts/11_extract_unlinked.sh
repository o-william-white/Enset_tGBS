#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_extract_unlinked

module load anaconda3
source activate stacks_workflow
ml python/3.6.3
ml R/3.6.1

mkdir -p extract_unlinked

# extract unlinked
# minimum difference between 0.0 and 1.0 to keep a second snp = 0.5 (recommended)
# maximum distance in base pairs to consider linked SNPs = 10kb 
python \
  additional_scripts/stacks_workflow/00-scripts/11_extract_unlinked_snps_genome.py \
  snp_duplication/populations.snps.singleton.vcf \
  0.5 \
  10000 \
  extract_unlinked/populations.snps.unlinked.vcf

# count the number of singleton snps
N_SING=$(grep -e "^#" -v -c snp_duplication/populations.snps.singleton.vcf)

# count the number of snps retained
N_SNPS=$(grep -e "^#" -v -c extract_unlinked/populations.snps.unlinked.vcf)

# print the number of snps
echo ${N_SING} singleton snps
echo ${N_SNPS} unlinked snps retained

echo done

