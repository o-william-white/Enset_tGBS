#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_bwa_index

module add bwa

# cp ref to dir
cp /data/scratch/mpx469/tGBS_enset_project/bedadeti_assembly/GCA_000818735.3_Bedadeti_annotated/GCA_000818735.3_Bedadeti_annotated_genomic.fna.gz .

# gunzip
gunzip GCA_000818735.3_Bedadeti_annotated_genomic.fna.gz

# index the reference genome
bwa index GCA_000818735.3_Bedadeti_annotated_genomic.fna

echo done

