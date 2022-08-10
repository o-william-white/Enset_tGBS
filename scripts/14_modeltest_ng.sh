#!/bin/bash
#$ -pe smp 8
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_modeltest_ng

export PATH=/data/home/mpx469/software/modeltest/modeltest-ng-0.1.5/build/bin/:$PATH

# mk output dir
mkdir -p modeltest_ng

# modeltest_ng
modeltest-ng \
   --datatype nt \
   --input convert_alignment_format/populations.all.seq.phylip \
   --output modeltest_ng/modeltest_ng \
   --processes 4

