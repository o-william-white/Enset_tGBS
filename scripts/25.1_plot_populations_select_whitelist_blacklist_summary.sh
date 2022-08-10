
ml R/3.6.1

# get summary stats for each assembly
echo -e "length\tloci\tsites\tfiltered\tvariant" > populations/summary_populations_select_whitelist_blacklist.txt

for LENGTH in `seq 70 10 120`; do 
    echo -e "${LENGTH}\t$(grep Kept populations/populations_select_${LENGTH}_all_snps_whitelist_blacklist/populations.log | cut -f 2,6,8,14 -d " " | sed -e 's/ /\t/g')" >> populations/summary_populations_select_whitelist_blacklist.txt
done

# summary plot
Rscript additional_scripts/plot_populations_summary.R populations/summary_populations_select_whitelist_blacklist.txt

