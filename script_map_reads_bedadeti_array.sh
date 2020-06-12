#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_map_reads_bedadeti_array
#$ -t 1-1698
#$ -tc 100


module add bwa
module add samtools

LENGTH=$(sed -n "${SGE_TASK_ID}p" input_args | cut -f 1 -d " ")
SAMPLE=$(sed -n "${SGE_TASK_ID}p" input_args | cut -f 2 -d " ")

# cp process_radtags output to dir
cp /data/scratch/mpx469/tGBS_enset_project/process_radtags/process_radtags_${LENGTH}_output/${SAMPLE}.cutadapt.with.cutsites.fq.gz tmp.${SAMPLE}.${LENGTH}.fq.gz

# gunzip
gunzip tmp.${SAMPLE}.${LENGTH}.fq.gz

# bwa mem
bwa mem GCA_000818735.3_Bedadeti_annotated_genomic.fna tmp.${SAMPLE}.${LENGTH}.fq | samtools view -o bwa_mem_${LENGTH}_output/${SAMPLE}.bam - 

# rm tmp file
rm tmp.${SAMPLE}.${LENGTH}.fq

# filter out reads that have not mapped uniquely and sort
# XA - alternative hits
# SA - chimeric read
samtools view -h /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti/bwa_mem_${LENGTH}_output/${SAMPLE}.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools sort -o samtools_${LENGTH}_output/${SAMPLE}.unique.sorted.bam -T ${SAMPLE}.${LENGTH}.tmp -

# test if all mapped reads are unique
if [ "`samtools view samtools_${LENGTH}_output/${SAMPLE}.unique.sorted.bam |  wc -l`" == "`samtools view samtools_${LENGTH}_output/${SAMPLE}.unique.sorted.bam | cut -f 1 | sort | uniq | wc -l`" ]; then
   echo "all mapped reads are unique";
else
   echo "check read mapping!"
fi

# create index to allow viewing in igv
samtools index samtools_${LENGTH}_output/${SAMPLE}.unique.sorted.bam

echo done

