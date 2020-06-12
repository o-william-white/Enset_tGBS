#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=12:0:0
#$ -cwd
#$ -j y
#$ -N job_ftp_refseq_bacteria

# wget
wget -i ftp_refseq_bacteria/assembly_summary_select_paths -P ftp_refseq_bacteria/

# gunzip
gunzip ftp_refseq_bacteria/*.gz

# cat assemblies to single fasta
cat ftp_refseq_bacteria/*.fna > ftp_refseq_bacteria/refseq_bacteria.fasta

# check number of assemblies downloaded
echo number of assemblies downloaded = `ls -1 ftp_refseq_bacteria/*fna | wc -l`

# remove individual .fna files
rm ftp_refseq_bacteria/*.fna

echo done

