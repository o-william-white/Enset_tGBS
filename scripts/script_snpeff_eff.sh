#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_snpeff_eff

# snpeff
java -Xmx10g -jar snpEff.jar -classic bedadeti populations.snps.vcf > populations.snps.eff.vcf
 
# snpsift - wide
java -jar SnpSift.jar extractFields -s "," -e "." populations.snps.eff.vcf ID ALT "ANN[*].EFFECT" > snpsift_wide_eff.txt

# snpsift - long
cat populations.snps.eff.vcf | ./vcfEffOnePerLine.pl | grep -e "^##" -v | java -jar SnpSift.jar extractFields -s "," -e "." - ID ALT "ANN[*].EFFECT" > snpsift_long_eff.txt

