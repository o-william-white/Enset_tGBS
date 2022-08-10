#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_convert_alignment_format

module load python
source /data/home/mpx469/software/python-virtualenv/biopython-virtenv/bin/activate

# mkdir if not already present
mkdir -p convert_alignment_format

IN_DIR=populations/populations_80_all_snps_clean/
OUT_DIR=convert_alignment_format/

# remove line comment lines in stacks output
grep -e "^#" -v ${IN_DIR}/populations.all.phylip > ${OUT_DIR}/populations.all.phylip
grep -e "^#" -v ${IN_DIR}/populations.var.phylip > ${OUT_DIR}/populations.var.phylip

# python scripts
python additional_scripts/convert_alignment_format.py ${OUT_DIR}/ 

echo done

