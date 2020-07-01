
# run as:
# bash plot_plink_pca.sh 

for l in `seq 70 10 120`; do 
   for d in all_snps single_snp; do

      mkdir -p plink_${l}_${d}_blacklist_output

      plink --file /data/scratch/mpx469/tGBS_enset_project/populations/populations_${l}_${d}_blacklist_output/populations.plink --allow-extra-chr --maf 0.05 --pca --out plink_${l}_${d}_blacklist_output/plink

      Rscript plot_plink_pca.R plink_${l}_${d}_blacklist_output/

	done
done

