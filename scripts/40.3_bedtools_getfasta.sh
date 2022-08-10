#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_bedtools_getfasta

ml samtools
ml bedtools

mkdir -p bedtools

# cp ref to dir
cp bedadeti/GCA_000818735.3_Bedadeti_annotated_genomic.fna.gz bedtools/

# decompress
gunzip bedtools/GCA_000818735.3_Bedadeti_annotated_genomic.fna.gz

# create cp with differenct ending to keep as decompressed version 
cp bedtools/GCA_000818735.3_Bedadeti_annotated_genomic.fna bedtools/GCA_000818735.3_Bedadeti_annotated_genomic.fasta

# bgzip compress for samtools faidx compatability
bgzip bedtools/GCA_000818735.3_Bedadeti_annotated_genomic.fna

# samtools faidx to create genome file (.fai)
samtools faidx bedtools/GCA_000818735.3_Bedadeti_annotated_genomic.fna.gz

# cp bed files from bedtools
cp bayescan/outliers.bed bedtools/
cp bayescan/background.bed bedtools/

# add 10kb to either end
bedtools slop -i bedtools/outliers.bed   -g bedtools/GCA_000818735.3_Bedadeti_annotated_genomic.fna.gz.fai -b 10000 > bedtools/outliers_10kb.bed
bedtools slop -i bedtools/background.bed -g bedtools/GCA_000818735.3_Bedadeti_annotated_genomic.fna.gz.fai -b 10000 > bedtools/background_10kb.bed

# extract fasta files
bedtools getfasta -fi bedtools/GCA_000818735.3_Bedadeti_annotated_genomic.fasta -bed bedtools/outliers_10kb.bed   | fold -w 60 | tr [:lower:] [:upper:] > bedtools/outliers_10kb.fasta
bedtools getfasta -fi bedtools/GCA_000818735.3_Bedadeti_annotated_genomic.fasta -bed bedtools/background_10kb.bed | fold -w 60 | tr [:lower:] [:upper:] > bedtools/background_10kb.fasta

echo Complete!

