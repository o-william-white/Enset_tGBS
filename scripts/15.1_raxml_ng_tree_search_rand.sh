#!/bin/bash
#$ -pe smp 8
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -t 1-5000
#$ -tc 25
#$ -N job_raxml_ng_tree_search_rand

module add gcc/7.1.0
export PATH=/data/home/mpx469/software/raxml-ng-pthreads/raxml-ng/bin/:$PATH

SEED=`awk "NR==${SGE_TASK_ID}" raxml_ng/random_seed.txt`

raxml-ng \
   --msa raxml_ng/parse.raxml.rba \
   --search \
   --tree rand{1} \
   --threads 8 \
   --prefix raxml_ng/tree_search_rand_${SGE_TASK_ID} \
   --seed ${SEED}

echo done

