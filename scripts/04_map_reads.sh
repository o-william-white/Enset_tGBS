#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_map_reads_bedadeti
#$ -t 1-1698
#$ -tc 100

module add bwa
module add samtools

LENGTH=$(sed -n "${SGE_TASK_ID}p" input_args_process_radtags | cut -f 1 -d " ")
SAMPLE=$(sed -n "${SGE_TASK_ID}p" input_args_process_radtags | cut -f 2 -d " ")

# mkdir if not already present
mkdir -p map_reads_bedadeti
mkdir -p map_reads_bedadeti/bwa_mem_${LENGTH}
mkdir -p map_reads_bedadeti/samtools_${LENGTH}

# bwa mem
bwa mem \
  bedadeti/GCA_000818735.3_Bedadeti_annotated_genomic.fna.gz \
  process_radtags/process_radtags_${LENGTH}/${SAMPLE}.with.cutsites.fq.gz \
  | samtools view -o map_reads_bedadeti/bwa_mem_${LENGTH}/${SAMPLE}.bam - 

# filter out reads that have not mapped uniquely and sort
# XA - alternative hits
# SA - chimeric read
samtools view -h map_reads_bedadeti/bwa_mem_${LENGTH}/${SAMPLE}.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' \
  | samtools sort -o map_reads_bedadeti/samtools_${LENGTH}/${SAMPLE}.bam -T ${SAMPLE}.${LENGTH}.tmp -

# test if all mapped reads are unique
if [ "`samtools view map_reads_bedadeti/samtools_${LENGTH}/${SAMPLE}.bam |  wc -l`" == "`samtools view map_reads_bedadeti/samtools_${LENGTH}/${SAMPLE}.bam | cut -f 1 | sort | uniq | wc -l`" ]; then
   echo "all mapped reads are unique";
else
   echo "check read mapping!"
fi

# create index to allow viewing in igv
samtools index map_reads_bedadeti/samtools_${LENGTH}/${SAMPLE}.bam

echo done
