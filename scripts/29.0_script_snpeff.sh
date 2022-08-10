#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_snpeff
#$ -t 1-6

ml R/3.6.1

SNP_TYPE=$(echo -e "diverged\nduplicated\nhighcov\nlowconf\nmas\nsingleton" | sed -n "${SGE_TASK_ID}p")

# set dir
mkdir -p snpeff
mkdir -p snpeff/${SNP_TYPE}
mkdir -p snpeff/${SNP_TYPE}/data
mkdir -p snpeff/${SNP_TYPE}/data/genomes
mkdir -p snpeff/${SNP_TYPE}/data/bedadeti

# Get annotation files
wget -P snpeff/${SNP_TYPE}/data/bedadeti https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/818/735/GCA_000818735.3_Bedadeti_annotated/GCA_000818735.3_Bedadeti_annotated_genomic.gtf.gz
mv snpeff/${SNP_TYPE}/data/bedadeti/GCA_000818735.3_Bedadeti_annotated_genomic.gtf.gz snpeff/${SNP_TYPE}/data/bedadeti/genes.gtf.gz

# Get the genome
wget -P snpeff/${SNP_TYPE}/data/genomes https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/818/735/GCA_000818735.3_Bedadeti_annotated/GCA_000818735.3_Bedadeti_annotated_genomic.fna.gz
mv snpeff/${SNP_TYPE}/data/genomes/GCA_000818735.3_Bedadeti_annotated_genomic.fna.gz snpeff/${SNP_TYPE}/data/genomes/bedadeti.fa.gz

# cp scripts from software download to snpeff dir
cp /data/home/mpx469/software/SnpEff/snpEff/snpEff.jar snpeff/${SNP_TYPE}/
cp /data/home/mpx469/software/SnpEff/snpEff/SnpSift.jar snpeff/${SNP_TYPE}/
cp /data/home/mpx469/software/SnpEff/snpEff/snpEff.config snpeff/${SNP_TYPE}/
cp /data/home/mpx469/software/SnpEff/snpEff/scripts/vcfEffOnePerLine.pl snpeff/${SNP_TYPE}/

# add entry to config file
echo "bedadeti.genome : Enset" >> snpeff/${SNP_TYPE}/snpEff.config 

# ln populations input
ln -sf ${PWD}/populations/populations_select_80_all_snps_whitelist_blacklist_${SNP_TYPE}/populations.snps.vcf snpeff/${SNP_TYPE}/

# build database
cd snpeff/${SNP_TYPE}/
java -Xmx10g -jar snpEff.jar build -gtf22 -v bedadeti 

# snpeff
java -Xmx10g -jar snpEff.jar bedadeti populations.snps.vcf > populations.snps.ann.vcf
java -Xmx10g -jar snpEff.jar -classic bedadeti populations.snps.vcf > populations.snps.eff.vcf

# snpsift - wide
java -jar SnpSift.jar extractFields -s "," -e "." populations.snps.ann.vcf ID ALT "ANN[*].EFFECT" > snpsift_wide_ann.txt
java -jar SnpSift.jar extractFields -s "," -e "." populations.snps.eff.vcf ID ALT "ANN[*].EFFECT" > snpsift_wide_eff.txt

# snpsift - long
cat populations.snps.ann.vcf | ./vcfEffOnePerLine.pl | grep -e "^##" -v | java -jar SnpSift.jar extractFields -s "," -e "." - ID ALT "ANN[*].EFFECT" > snpsift_long_ann.txt
cat populations.snps.eff.vcf | ./vcfEffOnePerLine.pl | grep -e "^##" -v | java -jar SnpSift.jar extractFields -s "," -e "." - ID ALT "ANN[*].EFFECT" > snpsift_long_eff.txt

# format data
Rscript ../../additional_scripts/snpeff_format_annotation_table_ann.R
Rscript ../../additional_scripts/snpeff_format_annotation_table_eff.R

echo Complete!

