#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_plot_snp_summary_stats_blacklist
#$ -t 70,80,90,100,110,120
#$ -tc 6

ml R/3.6.1

mkdir -p snp_summary_plots

Rscript additional_scripts/plot_snp_summary_stats.R \
  populations/populations_${SGE_TASK_ID}_all_snps_blacklist/populations.snps.vcf \
  tGBS_metadata_phylogenetic_analysis.csv \
  snp_summary_plots/summary_all_snps_blacklist_${SGE_TASK_ID}

