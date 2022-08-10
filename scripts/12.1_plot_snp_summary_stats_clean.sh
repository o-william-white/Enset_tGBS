#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_plot_snp_summary_stats_clean

ml R/3.6.1

mkdir -p snp_summary_plots

Rscript additional_scripts/plot_snp_summary_stats.R \
  populations/populations_80_all_snps_clean/populations.snps.vcf \
  tGBS_metadata_phylogenetic_analysis.csv \
  snp_summary_plots/summary_all_snps_clean_80

