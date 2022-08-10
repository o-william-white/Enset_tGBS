#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=3G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_population_genetics

ml R/3.6.1

mkdir -p populations/populations_select_80_all_snps_whitelist_blacklist_concatenated 
ln -sf $PWD/concatenate_vcfs/populations.snps.vcf populations/populations_select_80_all_snps_whitelist_blacklist_concatenated

Rscript additional_scripts/pg_00_prop_het.R concatenated
Rscript additional_scripts/pg_01_maf.R concatenated
Rscript additional_scripts/pg_02_fis.R concatenated
Rscript additional_scripts/pg_03_allele_categories.R concatenated

echo done

