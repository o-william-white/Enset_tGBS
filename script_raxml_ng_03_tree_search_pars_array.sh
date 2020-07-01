#!/bin/bash
#$ -pe smp 8
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -t 1-5000
#$ -tc 25
#$ -N job_raxml_ng_03_tree_search_pars

module add gcc/7.1.0
export PATH=/data/home/mpx469/software/raxml-ng-pthreads/raxml-ng/bin/:$PATH

SEED=`awk "NR==${SGE_TASK_ID}" random_seed.txt`

raxml-ng \
   --msa raxml_ng_output/all_parse.raxml.rba \
   --search \
   --tree pars{1} \
   --threads 8 \
   --prefix raxml_ng_output/all_tree_search_pars_${SGE_TASK_ID} \
   --seed ${SEED} 

echo done

