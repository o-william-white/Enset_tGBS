#!/bin/bash
#$ -pe smp 16
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_iqtree
#$ -t 1-500
#$ -tc 50

export PATH=/data/home/mpx469/software/iq-tree/iqtree-1.6.12-Linux/bin/:$PATH

INPUT_ARGS=$(sed -n "${SGE_TASK_ID}p" input_args_iqtree)

mkdir -p iqtree

iqtree -s convert_alignment_format/populations.all.seq.phylip \
   -m MFP -bb 1000 -alrt 1000 -nt 16 -pre iqtree/iqtree ${INPUT_ARGS}

echo done

