#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G
#$ -j y
#$ -N job_missingness

module load vcftools
module load R

for l in `seq 70 10 120`; do 
   for d in all_snps single_snp; do
	   
      mkdir -p missingness_${l}_${d}_output

      vcftools --vcf /data/scratch/mpx469/tGBS_enset_project/populations/populations_${l}_${d}_blacklist_output/populations.snps.vcf \
         --missing-indv \
         --out missingness_${l}_${d}_output/missingness_${l}_${d}

      vcftools --vcf /data/scratch/mpx469/tGBS_enset_project/populations/populations_${l}_${d}_blacklist_output/populations.snps.vcf \
         --missing-site \
         --out missingness_${l}_${d}_output/missingness_${l}_${d}
     	 
      Rscript plot_missingness.R ${l} ${d}
	 
   done
done

echo done

