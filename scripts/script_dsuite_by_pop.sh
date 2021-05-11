#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_dsuite_by_pop

module load vcftools
module load bcftools
module add R
export PATH=/data/home/mpx469/software/Dsuite/Build/:$PATH

# Dsuite Dtrios
Dsuite Dtrios /data/scratch/mpx469/tGBS_enset_project/populations/populations_80_single_snp_blacklist_output/populations.snps.vcf by_pop.txt

echo done

