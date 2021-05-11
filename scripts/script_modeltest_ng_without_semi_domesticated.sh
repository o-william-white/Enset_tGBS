#!/bin/bash
#$ -pe smp 8
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_modeltest_ng_without_semi_domesticated

export PATH=/data/home/mpx469/software/modeltest/modeltest-ng-0.1.5/build/bin/:$PATH

# mk output dir
mkdir -p modeltest_ng_80_single_snp_without_semi_domesticated

# modeltest_ng
modeltest-ng \
   --datatype nt \
   --input /data/scratch/mpx469/tGBS_enset_project/convert_file_formats/convert_file_formats_80_single_snp/populations.all.seq.without.semi.domesticated.phylip \
   --output modeltest_ng_80_single_snp_without_semi_domesticated/modeltest_ng \
   --processes 4

