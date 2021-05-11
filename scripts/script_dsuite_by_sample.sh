#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_dsuite_by_sample

module load vcftools
module load bcftools
module add R
export PATH=/data/home/mpx469/software/Dsuite/Build/:$PATH
module load python/3.6.3
source /data/home/mpx469/software/numpy/bin/activate

# Dsuite Dtrios
Dsuite Dtrios /data/scratch/mpx469/tGBS_enset_project/populations/populations_80_single_snp_blacklist_output/populations.snps.vcf by_sample.txt -t raxml.renamed.newick

# Dsuite Fbranch
Dsuite Fbranch raxml.renamed.newick by_sample__tree.txt > by_sample_Fbranch.txt

# plot function
python /data/home/mpx469/software/Dsuite/utils/dtools.py by_sample_Fbranch.txt raxml.renamed.newick

echo done

