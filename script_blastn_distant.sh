#!/bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l h_vmem=1G
#$ -j y
#$ -l h_rt=1:0:0
#$ -N job_blastn_distant
#$ -t 1-12
#$ -tc 12

module load blast+
module load R

LENGTH=$(sed -n "${SGE_TASK_ID}p" input_args | cut -f 1 -d " " )
DATASET=$(sed -n "${SGE_TASK_ID}p" input_args | cut -f 2 -d " " )

# blastn
blastn \
   -query /data/scratch/mpx469/tGBS_enset_project/populations/populations_distant_${LENGTH}_${DATASET}_output/populations.loci.fa \
   -db organelle_seq.db \
   -out blast_out_distant_${LENGTH}_${DATASET} \
   -evalue 1e-30 \
   -outfmt 6 \
   -num_threads 16

# get top hits
Rscript top_hit.R blast_out_distant_${LENGTH}_${DATASET}

echo done

