#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_convert_file_formats_array
#$ -t 1-12

LENGTH=$(sed -n "${SGE_TASK_ID}p" input_args | cut -f 1 -d " " )
DATASET=$(sed -n "${SGE_TASK_ID}p" input_args | cut -f 2 -d " " )

# set paths
module load python
source /data/home/mpx469/software/python-virtualenv/biopython-virtenv/bin/activate

# create output directory
mkdir -p convert_file_formats_${LENGTH}_${DATASET}

INPUT=/data/scratch/mpx469/tGBS_enset_project/populations/populations_${LENGTH}_${DATASET}_blacklist_output/
OUTPUT=convert_file_formats_${LENGTH}_${DATASET}/

# remove line comment lines in stacks output
grep -e "^#" -v ${INPUT}/populations.all.phylip > ${OUTPUT}/populations.all.phylip
grep -e "^#" -v ${INPUT}/populations.var.phylip > ${OUTPUT}/populations.var.phylip

# python scripts
python convert_alignment_format.py ${OUTPUT}/ 

echo done

