#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_wget_and_index_bedadeti

module add bwa

# wget bedadeti genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/818/735/GCA_000818735.3_Bedadeti_annotated/GCA_000818735.3_* -P bedadeti/

# index
bwa index bedadeti/GCA_000818735.3_Bedadeti_annotated_genomic.fna.gz

echo done

