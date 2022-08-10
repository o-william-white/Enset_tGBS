#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_format_bayescan_input

module load vcftools
module load bcftools
module load samtools

mkdir -p bayescan

# get clone corrected domesticated and wild samples
# mlg_farthest_bitwise_monophyletic_single_rep.txt copied from local dir
grep -e 'Wild' -e 'Domesticated' estimate_mlg/mlg_farthest_bitwise_monophyletic_single_rep.txt  | cut -f 1 > bayescan/samples_domesticated_wild.txt

# ln vcf
# cp radpainter output
ln -sf ${PWD}/populations/populations_80_all_snps_clean/populations.snps.vcf ${PWD}/bayescan/

# sort vcf file based on chr and pos
# https://www.biostars.org/p/299659/
cat bayescan/populations.snps.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > bayescan/populations.snps.sorted.vcf

# bgzip compress and index with tabix
# https://www.biostars.org/p/59492/
bgzip -c bayescan/populations.snps.sorted.vcf > bayescan/populations.snps.sorted.vcf.gz
tabix -p vcf bayescan/populations.snps.sorted.vcf.gz

# select domesticated and wild samples for vcf file
bcftools view -Ov -S bayescan/samples_domesticated_wild.txt bayescan/populations.snps.sorted.vcf.gz > bayescan/populations.snps.sorted.dw.vcf

# maf filter of 0.05
vcftools --vcf bayescan/populations.snps.sorted.dw.vcf --maf 0.05 --recode --out bayescan/populations.snps.sorted.dw.maf

# create bayescan definitions file
grep -e 'Wild' -e 'Domesticated' estimate_mlg/mlg_farthest_bitwise_monophyletic_single_rep.txt | cut -f 1,4 > bayescan/vcf_to_bayescan_definitions.txt

# run pgdspider
cd bayescan
java -Xmx10g -Xms512M -jar /data/home/mpx469/software/pgd-spider/PGDSpider_2.1.1.5/PGDSpider2-cli.jar \
   -inputfile populations.snps.sorted.dw.maf.recode.vcf \
   -inputformat VCF \
   -outputfile populations.snps.maf.bayescan.txt \
   -outputrmat GESTE_BAYE_SCAN \
   -spid ../additional_scripts/vcf_to_bayescan.spid

echo Complete!

