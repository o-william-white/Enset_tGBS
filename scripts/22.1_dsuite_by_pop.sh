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
cd dsuite
Dsuite Dtrios populations.snps.vcf by_pop.txt

echo Complete!

