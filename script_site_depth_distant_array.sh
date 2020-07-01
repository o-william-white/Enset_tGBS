#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_site_depth_distant_array
#$ -t 1-336
#$ -tc 50

LENGTH=$(sed -n "${SGE_TASK_ID}p" input_args_site_depth_distant | cut -f 1 -d " " )
DATASET=$(sed -n "${SGE_TASK_ID}p" input_args_site_depth_distant | cut -f 2 -d " " )
SAMPLE=$(sed -n "${SGE_TASK_ID}p" input_args_site_depth_distant | cut -f 3 -d " " )

module load bcftools
module load vcftools

mkdir -p site_depth_distant_${LENGTH}_${DATASET}_output

# subset for sample and site depth
bcftools view -Ov --samples ${SAMPLE} vcf_sorted_indexed/populations_distant_${LENGTH}_${DATASET}_sorted.vcf.gz | vcftools --vcf - --site-depth --out site_depth_distant_${LENGTH}_${DATASET}_output/${SAMPLE}

echo done

