
# run as
# bash bash_blacklist_summary_stats.sh

for l in `seq 70 10 120`; do 

   # set paths for files      
   BLAST=blastn/blast_${l}_all_snps/blast_blacklist.txt
   DUPS=duplicate_loci/duplicate_loci_${l}_all_snps_output/duplicate_blacklist.txt
   DEPTH=locus_depth/loci_depth_${l}_all_snps_output/site_depth_blacklist.txt

   # cat and get uniq values
   cat ${BLAST} ${DUPS} ${DEPTH} | sort | uniq > blacklist_${l}_all_snps

   # summary counts
   BLAST_COUNT=$(wc -l ${BLAST} | cut -f 1 -d " ")   
   DUPS_COUNT=$(wc -l ${DUPS}   | cut -f 1 -d " ")
   DEPTH_COUNT=$(wc -l ${DEPTH} | cut -f 1 -d " ")
   SHARED_COUNT=$(cat ${BLAST} ${DUPS} ${DEPTH} | sort | uniq -dc | wc -l)
   TOTAL_COUNT=$(cat ${BLAST} ${DUPS} ${DEPTH} | sort | uniq | wc -l)

   # print summary counts
   echo type count             >  summary_blacklist_${l}_all_snps
   echo blast ${BLAST_COUNT}   >> summary_blacklist_${l}_all_snps
   echo dups ${DUPS_COUNT}     >> summary_blacklist_${l}_all_snps
   echo depth ${DEPTH_COUNT}   >> summary_blacklist_${l}_all_snps
   echo shared ${SHARED_COUNT} >> summary_blacklist_${l}_all_snps
   echo total ${TOTAL_COUNT}   >> summary_blacklist_${l}_all_snps
   
done

