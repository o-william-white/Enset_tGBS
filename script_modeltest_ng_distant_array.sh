#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_modeltest_ng_distant_array
#$ -t 70,80,90,100,110,120

export PATH=/data/home/mpx469/software/modeltest/modeltest-ng-0.1.5/build/bin/:$PATH

# mk output dir
mkdir -p modeltest_ng_distant_${SGE_TASK_ID}_single_snp

# modeltest_ng
modeltest-ng \
   --datatype nt \
   --input /data/scratch/mpx469/tGBS_enset_project/convert_file_formats/convert_file_formats_distant_${SGE_TASK_ID}_single_snp/populations.all.seq.phylip \
   --output modeltest_ng_distant_${SGE_TASK_ID}_single_snp/modeltest_ng \
   --processes 4

