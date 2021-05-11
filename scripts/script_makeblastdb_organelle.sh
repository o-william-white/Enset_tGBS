#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_makeblastdb

module load blast+

# makebdb
makeblastdb -in organelle.fasta -out organelle.db -dbtype nucl

echo done

