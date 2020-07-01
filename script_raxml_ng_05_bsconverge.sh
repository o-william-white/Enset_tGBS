#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_all_bsconverge

module add gcc/7.1.0
export PATH=/data/home/mpx469/software/raxml-ng-pthreads/raxml-ng/bin/:$PATH

# cat all bootstrap trees
cat raxml_ng_output/all_bootstrap_*.raxml.bootstraps > raxml_ng_output/all_bootstraps

# bscoverge
raxml-ng \
   --bsconverge \
   --bs-trees raxml_ng_output/all_bootstraps \
   --prefix raxml_ng_output/all_bsconvergence \
   --bs-cutoff 0.03 \
   --threads 4

