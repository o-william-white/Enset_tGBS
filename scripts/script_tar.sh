#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=24:0:0
#$ -cwd
#$ -j y
#$ -N job_tar

tar xvzf Data2Bio_final.tar.gz 

echo done

