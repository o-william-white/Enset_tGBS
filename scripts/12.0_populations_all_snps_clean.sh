#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=48G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_populations_all_snps_clean

module add stacks/2.41

# mkdir if not already present
mkdir -p populations
mkdir -p populations/populations_80_all_snps_clean

# get locus, chromosome, bp and col positions from original vcf
grep -e "^#" -v populations/populations_80_all_snps_blacklist/populations.sumstats.tsv | cut -f 1-4 | sort -k 1 -n | uniq > populations/populations_80_all_snps_clean/reference_locus_chr_bp_col.txt

# get chromosome and bp of clean snps
grep -e "^#" -v extract_unlinked/populations.snps.unlinked.vcf | cut -f 1,2 > populations/populations_80_all_snps_clean/clean_chr_bp.txt

# get locus and col of clean snps from reference data
# note that I have to increase the memory requirement for this step
grep -w -f populations/populations_80_all_snps_clean/clean_chr_bp.txt populations/populations_80_all_snps_clean/reference_locus_chr_bp_col.txt | cut -f 1,4 > populations/populations_80_all_snps_clean/whitelist.txt

echo Number of sites unlinked vcf: $(grep -e "^#" -v -c extract_unlinked/populations.snps.unlinked.vcf )
echo Number of sites in whitelist: $(wc -l populations/populations_80_all_snps_clean/whitelist.txt | cut -f 1 -d " ")

# rm intermediate files
rm populations/populations_80_all_snps_clean/reference_locus_chr_bp_col.txt
rm populations/populations_80_all_snps_clean/clean_chr_bp.txt

# populations
populations \
   -P gstacks/gstacks_80/ \
   -O populations/populations_80_all_snps_clean \
   -M gstacks/popmap.txt \
   -W populations/populations_80_all_snps_clean/whitelist.txt \
   --vcf --phylip --phylip-var --phylip-var-all --fasta-loci --plink --radpainter --ordered-export

echo Number of sites unlinked vcf: $(grep -e "^#" -v -c extract_unlinked/populations.snps.unlinked.vcf )
echo Number of sites clean vcf: $(grep -e "^#" -v -c populations/populations_80_all_snps_clean/populations.snps.vcf )

echo done

