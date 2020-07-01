#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_raxml_ng_1_parse

module add gcc/7.1.0
export PATH=/data/home/mpx469/software/raxml-ng-pthreads/raxml-ng/bin/:$PATH

# mk output dir
mkdir -p raxml_ng_output/

# grep log file to find the selected model
MODEL=$(grep Summary -A 7 /data/scratch/mpx469/tGBS_enset_project/modeltest_ng/modeltest_ng_80_single_snp/modeltest_ng.log | grep AICc | awk ' { print $2 } ')

raxml-ng \
   --msa /data/scratch/mpx469/tGBS_enset_project/convert_file_formats/convert_file_formats_80_single_snp/populations.all.seq.phylip \
   --model ${MODEL} \
   --parse \
   --prefix raxml_ng_output/all_parse

# generate random seed 
shuf -i 1-10000 -n 5000 > random_seed.txt

