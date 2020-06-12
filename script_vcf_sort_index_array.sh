#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_vcf_sort_index_array
#$ -t 1-12

module load bcftools
module load vcftools

LENGTH=$(sed -n "${SGE_TASK_ID}p" input_args | cut -f 1 -d " ")
DATASET=$(sed -n "${SGE_TASK_ID}p" input_args | cut -f 2 -d " ")

# cp vcf file  
cp /data/scratch/mpx469/tGBS_enset_project/populations/populations_${LENGTH}_${DATASET}_output/populations.snps.vcf vcf_sorted_indexed/populations_${LENGTH}_${DATASET}.vcf

# remove .unique.sorted endings from sample names
sed -i 's/.unique.sorted//g' vcf_sorted_indexed/populations_${LENGTH}_${DATASET}.vcf

# sort vcf file based on chr and pos
# https://www.biostars.org/p/299659/
cat vcf_sorted_indexed/populations_${LENGTH}_${DATASET}.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > vcf_sorted_indexed/populations_${LENGTH}_${DATASET}_sorted.vcf

# bgzip compress and index with tabix
# https://www.biostars.org/p/59492/
bgzip -c vcf_sorted_indexed/populations_${LENGTH}_${DATASET}_sorted.vcf > vcf_sorted_indexed/populations_${LENGTH}_${DATASET}_sorted.vcf.gz
tabix -p vcf vcf_sorted_indexed/populations_${LENGTH}_${DATASET}_sorted.vcf.gz

# rm unnecessary files
rm vcf_sorted_indexed/populations_${LENGTH}_${DATASET}_sorted.vcf
rm vcf_sorted_indexed/populations_${LENGTH}_${DATASET}.vcf

echo done

