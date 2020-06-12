
# run as
# bash script_06_write_summary_tables.sh

# paste columns together in specified order
# .tmp suffix used as intermeidate to add column names
for i in `seq 70 10 120`; do 
   paste /data/scratch/mpx469/tGBS_enset_project/sample_list.txt \
         count_00_raw.txt \
         count_01_trimmomatic.txt \
         count_02_cutadapt.txt \
         count_03_process_radtags_${i}.txt \
         count_04_bwa_mem_${i}_mapped.txt \
         count_04_bwa_mem_${i}_unmapped.txt \
         count_04_bwa_mem_${i}_XA_SA.txt \
         count_05_samtools_${i}_mapped.txt \
         count_05_samtools_${i}_unmapped.txt > summary_table_${i}.tmp
done

# create tab delimited column names
echo -e "sequence_id\traw\ttrimmomatic\tcutadapt\tprocess_radtags\tbwa_mem_mapped\tbwa_mem_unmapped\tbwa_mem_XA_SA\tsamtools_mapped\tsamtools_unmapped" > colnames

# cat colnames and summary .tmp summary table
for i in `seq 70 10 120`; do
   cat colnames summary_table_${i}.tmp > summary_table_${i}
done

# rm tmp summary file and colmanes
rm summary_table_*.tmp
rm colnames


