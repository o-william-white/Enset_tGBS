#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G
#$ -N job-bwa-index

module add bwa

# index the reference genome
bwa index GCA_000331365.3_Ensete_JungleSeeds_v3.0_genomic.fna

