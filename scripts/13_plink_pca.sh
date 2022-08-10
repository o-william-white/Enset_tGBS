
ml plink
ml R/3.6.1

mkdir -p plink_pca

# need a tab-delimited text file with family IDs in the first column and within-family IDs for outgroup samples
grep Outgroup tGBS_metadata_phylogenetic_analysis.csv | cut -f2 -d "," | grep -f - populations/populations_80_all_snps_clean/populations.plink.ped | cut -f 1,2 > plink_pca/outgroups.txt

# plink
plink \
  --file populations/populations_80_all_snps_clean/populations.plink \
  --allow-extra-chr \
  --maf 0.05 \
  --remove plink_pca/outgroups.txt \
  --pca \
  --out plink_pca/plink

# plot
Rscript additional_scripts/plot_plink_pca.R plink_pca

