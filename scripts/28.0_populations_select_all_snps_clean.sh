#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=250G
#$ -l highmem
#$ -l h_rt=12:0:0
#$ -cwd
#$ -j y
#$ -N job_populations_all_snps_whitelist_blacklist_snp_type
#$ -t 1-6

module add stacks/2.41

SNP_TYPE=$(echo -e "diverged\nduplicated\nhighcov\nlowconf\nmas\nsingleton" | sed -n "${SGE_TASK_ID}p")

# mkdir if not already present
mkdir -p populations/populations_select_80_all_snps_whitelist_blacklist_${SNP_TYPE}

# get locus, chromosome, bp and col positions from original vcf
grep -e "^#" -v populations/populations_select_80_all_snps_whitelist_blacklist/populations.sumstats.tsv | cut -f 1-4 | sort -k 1 -n | uniq > populations/populations_select_80_all_snps_whitelist_blacklist_${SNP_TYPE}/reference_locus_chr_bp_col.txt

# get chromosome and bp of target snp type
grep -e "^#" -v extract_unlinked_select/populations.snps.${SNP_TYPE}.unlinked.vcf | cut -f 1,2 > populations/populations_select_80_all_snps_whitelist_blacklist_${SNP_TYPE}/${SNP_TYPE}_chr_bp.txt

# get locus and col of snp type from reference data
# note that I have to increase the memory requirement for this step
grep -w -f populations/populations_select_80_all_snps_whitelist_blacklist_${SNP_TYPE}/${SNP_TYPE}_chr_bp.txt populations/populations_select_80_all_snps_whitelist_blacklist_${SNP_TYPE}/reference_locus_chr_bp_col.txt | cut -f 1,4 > populations/populations_select_80_all_snps_whitelist_blacklist_${SNP_TYPE}/whitelist.txt

echo Number of sites unlinked vcf: $(grep -e "^#" -v -c extract_unlinked_select/populations.snps.${SNP_TYPE}.unlinked.vcf )
echo Number of sites in whitelist: $(wc -l populations/populations_select_80_all_snps_whitelist_blacklist_${SNP_TYPE}/whitelist.txt | cut -f 1 -d " ")

# rm intermediate files
rm populations/populations_select_80_all_snps_whitelist_blacklist_${SNP_TYPE}/reference_locus_chr_bp_col.txt
rm populations/populations_select_80_all_snps_whitelist_blacklist_${SNP_TYPE}/${SNP_TYPE}_chr_bp.txt

# populations
populations \
   -P gstacks/gstacks_80/ \
   -O populations/populations_select_80_all_snps_whitelist_blacklist_${SNP_TYPE} \
   -M gstacks/popmap_select_all.txt \
   -W populations/populations_select_80_all_snps_whitelist_blacklist_${SNP_TYPE}/whitelist.txt \
   --vcf --phylip --phylip-var --phylip-var-all --fasta-loci --plink --radpainter --ordered-export

echo Number of sites unlinked vcf: $(grep -e "^#" -v -c extract_unlinked_select/populations.snps.${SNP_TYPE}.unlinked.vcf )
echo Number of sites ${SNP_TYPE} vcf: $(grep -e "^#" -v -c populations/populations_select_80_all_snps_whitelist_blacklist_${SNP_TYPE}/populations.snps.vcf )

echo done

