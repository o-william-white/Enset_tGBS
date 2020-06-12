#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=30G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_makeblastdb

module load blast+

# makebdb
makeblastdb -in ftp_refseq_bacteria/refseq_bacteria.fasta -out ftp_refseq_bacteria/refseq_bacteria.db -dbtype nucl 

#-parse_seqids -taxid_map ftp_refseq_bacteria/taxa_map.txt

echo done

