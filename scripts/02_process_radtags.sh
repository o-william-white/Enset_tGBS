#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_process_radtags
#$ -t 1-1698
#$ -tc 100

module load use.dev
module load stacks/2.41

LENGTH=$(sed -n "${SGE_TASK_ID}p" input_args_process_radtags | cut -f 1 -d " ")
SAMPLE=$(sed -n "${SGE_TASK_ID}p" input_args_process_radtags | cut -f 2 -d " ")

# mkdir if not already present
mkdir -p process_radtags
mkdir -p process_radtags_${LENGTH}

process_radtags \
   --disable-rad-check \
   -f cutadapt/${SAMPLE}.with.cutsites.fq.gz \
   -t ${LENGTH} \
   --len-limit ${LENGTH} \
   -o process_radtags/process_radtags_${LENGTH}/

echo done

