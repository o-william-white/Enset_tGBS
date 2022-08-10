#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_iq_tree_write_args

# write input args
for i in `seq 0.1 0.1 0.5`; do
   for r in {1..100}; 
      do echo -pers ${i} -nstop 1000 -pre iqtree/out_pers_${i}_r${r}; 
   done 
done > input_args_iqtree

echo done

