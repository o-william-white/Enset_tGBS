#!/bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l h_vmem=1G
#$ -j y
#$ -l h_rt=1:0:0
#$ -N job_blastn_select
#$ -t 70,80,90,100,110,120

ml blast+
ml R/3.6.1

# mkdir if not already present
mkdir -p blacklist/select_blastn/
mkdir -p blacklist/select_blastn/blastn_${SGE_TASK_ID}

# blastn
blastn \
   -query populations/populations_select_${SGE_TASK_ID}_all_snps_whitelist/populations.loci.fa \
   -db blacklist/blastn/organelle.db \
   -out blacklist/select_blastn/blastn_${SGE_TASK_ID}/blast_out \
   -evalue 1e-30 \
   -outfmt 6 \
   -num_threads 16

# get top hits
Rscript additional_scripts/top_hit.R blacklist/select_blastn/blastn_${SGE_TASK_ID}/blast_out

# write blacklists
cut -f 1 blacklist/select_blastn/blastn_${SGE_TASK_ID}/blast_out_top_hits | sed 's/CLocus_//g' | grep qseqid -v > blacklist/select_blastn/blastn_${SGE_TASK_ID}/blacklist.txt

echo done

