#!/bin/bash
#$ -pe smp 16
#$ -l h_vmem=1G
#$ -l h_rt=12:0:0
#$ -cwd
#$ -j y
#$ -N job_iq_tree
#$ -t 1-1000
#$ -tc 50

export PATH=/data/home/mpx469/software/iq-tree/iqtree-1.6.12-Linux/bin/:$PATH

INPUT_ARGS=$(sed -n "${SGE_TASK_ID}p" input_args)

mkdir -p iq_tree_output

iqtree -s /data/scratch/mpx469/tGBS_enset_project/convert_file_formats/convert_file_formats_80_single_snp/populations.all.seq.phylip \
   -m MFP -bb 1000 -alrt 1000 -nt 16 ${INPUT_ARGS}

echo done

