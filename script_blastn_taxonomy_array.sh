#!/bin/sh
#$ -cwd
#$ -j y
#$ -pe smp 16
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -t 1-283
#$ -tc 50
#$ -N job_blastn_taxonomy_array

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" /data/scratch/mpx469/tGBS_enset_project/sample_list.txt)

module load seqtk
export PATH=/data/home/mpx469/software/cdhit/cd-hit-auxtools/:$PATH
module load blast+
module load R

# mk output dir if not already present
mkdir -p seqtk_output
mkdir -p cd_hit_output
mkdir -p blastn_output
mkdir -p summary_table_output/
mkdir -p summary_table_output/species/
mkdir -p summary_table_output/genera/

# convert trimmomatic output to fasta
seqtk seq -a /data/scratch/mpx469/tGBS_enset_project/trimmomatic/trimmomatic_output/${INPUT_FILE}.trimmomatic.fq.gz > seqtk_output/${INPUT_FILE}.fa

# collapse duplicate sequences
cd-hit-dup -i seqtk_output/${INPUT_FILE}.fa -o cd_hit_output/${INPUT_FILE}_cd_hit.fa 

# blastn
blastn \
   -db /data/SBCS-Ethiopia/tGBS_enset_project/refseq_bacteria/ftp_refseq_bacteria/refseq_bacteria_v2.db \
   -query cd_hit_output/${INPUT_FILE}_cd_hit.fa \
   -outfmt "6 std staxid ssciname" \
   -max_target_seqs 10 \
   -max_hsps 1 \
   -evalue 1e-25 \
   -out blastn_output/${INPUT_FILE}_blastn_out \
   -num_threads 16

# identify taxonomy
Rscript get_taxonomy.R ${INPUT_FILE}

echo done

