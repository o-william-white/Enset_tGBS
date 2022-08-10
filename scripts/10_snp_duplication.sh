#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_snp_duplication

module load anaconda3
source activate stacks_workflow
#ml python/3.6.3
ml R/3.6.1

mkdir -p snp_duplication

# rm link if already present
rm -f snp_duplication/populations.snps.vcf

# symbolic link
ln -s /data/scratch/mpx469/tGBS_enset_project/populations/populations_80_all_snps_blacklist/populations.snps.vcf snp_duplication/

# get snp info
python3 \
  additional_scripts/stacks_workflow/00-scripts/08_extract_snp_duplication_info.py \
  snp_duplication/populations.snps.vcf \
  snp_duplication/populations.snps.all.stats

# plot
Rscript additional_scripts/stacks_workflow/00-scripts/09_classify_snps.R \
  snp_duplication/populations.snps.all.stats

# split vcf
python3 \
  additional_scripts/stacks_workflow/00-scripts/10_split_vcf_in_categories.py \
  snp_duplication/populations.snps.vcf \
  snp_duplication/populations.snps.all.stats.categorized

echo done

