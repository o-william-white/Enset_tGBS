#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_plot_snp_summary_stats_select_whitelist_blacklist
#$ -t 70,80,90,100,110,120
#$ -tc 6

ml R/3.6.1

mkdir -p snp_summary_plots_select

Rscript additional_scripts/plot_snp_summary_stats.R \
  populations/populations_select_${SGE_TASK_ID}_all_snps_whitelist_blacklist/populations.snps.vcf \
  tGBS_metadata_population_genetics.csv \
  snp_summary_plots_select/summary_select_all_snps_whitelist_blacklist_${SGE_TASK_ID}

echo Complete!

