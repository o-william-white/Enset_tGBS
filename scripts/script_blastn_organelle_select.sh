#!/bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l h_vmem=1G
#$ -j y
#$ -l h_rt=1:0:0
#$ -N job_blastn_select
#$ -t 1-12
#$ -tc 12

module load blast+
module load R

LENGTH=$(sed -n "${SGE_TASK_ID}p" input_args | cut -f 1 -d " " )
DATASET=$(sed -n "${SGE_TASK_ID}p" input_args | cut -f 2 -d " " )

# mkdir if not already present
mkdir -p blastn_select_${LENGTH}_${DATASET}

# blastn
blastn \
   -query /data/scratch/mpx469/tGBS_enset_project/populations/populations_select_${LENGTH}_${DATASET}_whitelist_output/populations.loci.fa \
   -db organelle.db \
   -out blastn_select_${LENGTH}_${DATASET}/blast_out \
   -evalue 1e-30 \
   -outfmt 6 \
   -num_threads 16

# get top hits
Rscript top_hit.R blastn_select_${LENGTH}_${DATASET}/blast_out

# write blacklists
cut -f 1 blastn_select_${LENGTH}_${DATASET}/blast_out_top_hits | sed 's/CLocus_//g' | grep qseqid -v > blastn_select_${LENGTH}_${DATASET}/blastn_blacklist.txt

echo done

