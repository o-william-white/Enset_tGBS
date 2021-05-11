#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_format_bayescan_input

module load vcftools

# cp vcf to dir
cp /data/scratch/mpx469/tGBS_enset_project/populations/populations_select_80_all_snps_whitelist_blacklist_output/populations.snps.vcf .

# maf filter of 0.05
vcftools --vcf populations.snps.vcf --maf 0.05 --recode --out populations.snps.maf.vcf

# run pgdspider
java -Xmx10g -Xms512M -jar /data/home/mpx469/software/pgd-spider/PGDSpider_2.1.1.5/PGDSpider2-cli.jar \
   -inputfile populations.snps.maf.vcf.recode.vcf \
   -inputformat VCF \
   -outputfile populations.snps.maf.bayescan.txt \
   -outputrmat GESTE_BAYE_SCAN \
   -spid vcf_to_bayescan.spid

