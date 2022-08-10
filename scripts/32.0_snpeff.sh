#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_snpeff

ml R/3.6.1

SNP_TYPE=$(echo -e "diverged\nduplicated\nhighcov\nlowconf\nmas\nsingleton" | sed -n "${SGE_TASK_ID}p")

# set dir
mkdir -p snpeff
mkdir -p snpeff/concatenated
mkdir -p snpeff/concatenated/data
mkdir -p snpeff/concatenated/data/genomes
mkdir -p snpeff/concatenated/data/bedadeti

# Get annotation files
wget -P snpeff/concatenated/data/bedadeti https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/818/735/GCA_000818735.3_Bedadeti_annotated/GCA_000818735.3_Bedadeti_annotated_genomic.gtf.gz
mv snpeff/concatenated/data/bedadeti/GCA_000818735.3_Bedadeti_annotated_genomic.gtf.gz snpeff/concatenated/data/bedadeti/genes.gtf.gz

# Get the genome
wget -P snpeff/concatenated/data/genomes https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/818/735/GCA_000818735.3_Bedadeti_annotated/GCA_000818735.3_Bedadeti_annotated_genomic.fna.gz
mv snpeff/concatenated/data/genomes/GCA_000818735.3_Bedadeti_annotated_genomic.fna.gz snpeff/concatenated/data/genomes/bedadeti.fa.gz

# cp scripts from software download to snpeff dir
cp /data/home/mpx469/software/SnpEff/snpEff/snpEff.jar snpeff/concatenated/
cp /data/home/mpx469/software/SnpEff/snpEff/SnpSift.jar snpeff/concatenated/
cp /data/home/mpx469/software/SnpEff/snpEff/snpEff.config snpeff/concatenated/
cp /data/home/mpx469/software/SnpEff/snpEff/scripts/vcfEffOnePerLine.pl snpeff/concatenated/

# add entry to config file
echo "bedadeti.genome : Enset" >> snpeff/concatenated/snpEff.config 

# ln populations input
ln -sf ${PWD}/concatenate_vcfs/populations.snps.vcf snpeff/concatenated/

# build database
cd snpeff/concatenated/
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

