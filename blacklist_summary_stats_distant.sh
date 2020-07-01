
# run as
# bash blacklist_summary_stats_distant.sh

for l in `seq 70 10 120`; do 
   for d in all_snps single_snp; do

      BLAST=$(wc -l  blastn/blacklist_blast_distant_${l}_${d} | cut -f 1 -d " ")

      DUPS=$(wc -l duplicates/blacklist_duplicates_distant_${l}_${d} | cut -f 1 -d " ")

      DEPTH=$(wc -l site_depth/blacklist_site_depth_distant_${l}_${d} | cut -f 1 -d " ")

      SHARED=$(cat blastn/blacklist_blast_distant_${l}_${d} duplicates/blacklist_duplicates_distant_${l}_${d} site_depth/blacklist_site_depth_distant_${l}_${d} | sort | uniq -dc | wc -l)

      TOTAL=$(cat blastn/blacklist_blast_distant_${l}_${d} duplicates/blacklist_duplicates_distant_${l}_${d} site_depth/blacklist_site_depth_distant_${l}_${d} | sort | uniq | wc -l)

      echo blast ${BLAST}   >> summary_distant_${l}_${d}
      echo dups ${DUPS}     >> summary_distant_${l}_${d}
      echo depth ${DEPTH}   >> summary_distant_${l}_${d}
      echo shared ${SHARED} >> summary_distant_${l}_${d}
      echo total ${TOTAL}   >> summary_distant_${l}_${d}
   
      cat blastn/blacklist_blast_distant_${l}_${d} duplicates/blacklist_duplicates_distant_${l}_${d} site_depth/blacklist_site_depth_distant_${l}_${d} | sort | uniq > blacklist_overall_distant_${l}_${d}
   
   done
done

