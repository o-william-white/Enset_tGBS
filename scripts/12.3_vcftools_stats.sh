ml vcftools/

mkdir -p vcftools_stats

vcftools --vcf populations/populations_80_all_snps_clean/populations.snps.vcf --het --out vcftools_stats/output

