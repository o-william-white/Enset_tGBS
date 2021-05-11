#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_build_database

java -Xmx10g -jar snpEff.jar build -gtf22 -v bedadeti 

