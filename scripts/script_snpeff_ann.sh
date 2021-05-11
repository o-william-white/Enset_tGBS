#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_snpeff_ann

# snpeff
java -Xmx10g -jar snpEff.jar bedadeti populations.snps.vcf > populations.snps.ann.vcf

# snpsift - wide
java -jar SnpSift.jar extractFields -s "," -e "." populations.snps.ann.vcf ID ALT "ANN[*].EFFECT" > snpsift_wide_ann.txt

# snpsift - long
cat populations.snps.ann.vcf | ./vcfEffOnePerLine.pl | grep -e "^##" -v | java -jar SnpSift.jar extractFields -s "," -e "." - ID ALT "ANN[*].EFFECT" > snpsift_long_ann.txt
 
