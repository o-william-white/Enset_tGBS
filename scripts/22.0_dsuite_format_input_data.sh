
ml R/3.6.1

mkdir -p dsuite

# write sets
Rscript additional_scripts/dsuite_write_sets.R

# rename tip tables
Rscript additional_scripts/dsuite_write_tree_renamed_tip_labels.R

# ln populations output
ln -sf ${PWD}/populations/populations_80_all_snps_clean/populations.snps.vcf dsuite


