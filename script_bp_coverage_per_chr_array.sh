#!/bin/sh
#$ -cwd
#$ -j y
#$ -pe smp 16
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -t 1-283
#$ -tc 50
#$ -N job_bp_coverage_per_chr_array

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" /data/scratch/mpx469/tGBS_enset_project/sample_list.txt)

module load bedtools
module load samtools
module load bedops
module load R

# mk output dir if not already present
mkdir -p summary_table_output/xcm_cov
mkdir -p summary_table_output/xcm_bed

# run r script
Rscript bp_coverage_per_chr.R ${INPUT_FILE}

echo done

