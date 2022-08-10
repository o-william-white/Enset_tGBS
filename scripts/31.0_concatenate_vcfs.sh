#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_concatenate_vcfs

ml vcftools

mkdir -p concatenate_vcfs

# grep header from files to add to the singleton dataset
grep -e "^#" -v populations/populations_select_80_all_snps_whitelist_blacklist_lowconf/populations.snps.vcf > concatenate_vcfs/tmp_lowconf.vcf
grep -e "^#" -v populations/populations_select_80_all_snps_whitelist_blacklist_mas/populations.snps.vcf > concatenate_vcfs/tmp_mas.vcf

# concatenate singleton, lowconf and mas
cat populations/populations_select_80_all_snps_whitelist_blacklist_singleton/populations.snps.vcf \
  concatenate_vcfs/tmp_lowconf.vcf \
  concatenate_vcfs/tmp_mas.vcf | vcf-sort > concatenate_vcfs/populations.snps.vcf

# rm intermediate
rm concatenate_vcfs/tmp_*

echo Complete!

