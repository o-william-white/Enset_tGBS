
ml R/3.6.1

# get summary stats for loci assembled at each read length
echo -e "length\tloci\treads\tcoverage" > gstacks/summary_gstacks.txt

for LENGTH in `seq 70 10 120`; do 
  LOCI=$(grep Built gstacks/gstacks_${LENGTH}/gstacks.log | cut -f 2 -d " ")
  READS=$(grep Built gstacks/gstacks_${LENGTH}/gstacks.log | cut -f 5 -d " ") 
  COVERAGE=$(grep coverage gstacks/gstacks_${LENGTH}/gstacks.log | cut -f 6 -d " " | sed -e 's/mean=//g' -e 's/x,//g' )
  echo -e "${LENGTH}\t${LOCI}\t${READS}\t${COVERAGE}" >> gstacks/summary_gstacks.txt
done

# plot summary stats
Rscript additional_scripts/plot_gstacks_summary.R

