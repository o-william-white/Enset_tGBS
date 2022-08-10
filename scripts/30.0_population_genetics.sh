#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=3G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_population_genetics
#$ -t 1-6

ml R/3.6.1

SNP_TYPE=$(echo -e "diverged\nduplicated\nhighcov\nlowconf\nmas\nsingleton" | sed -n "${SGE_TASK_ID}p")

Rscript additional_scripts/pg_00_prop_het.R ${SNP_TYPE}
Rscript additional_scripts/pg_01_maf.R ${SNP_TYPE}
Rscript additional_scripts/pg_02_fis.R ${SNP_TYPE}
Rscript additional_scripts/pg_03_allele_categories.R ${SNP_TYPE}

echo done

