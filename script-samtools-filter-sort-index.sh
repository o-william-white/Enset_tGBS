#!/bin/sh
#$ -cwd
#$ -j y
#$ -l h_rt=2:0:0
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -N job-samtools
#$ -t 1-283
#$ -tc 36

module add samtools

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" /data/scratch/mpx469/sample-list.txt)

# filter out reads that have not mapped uniquely
# XA - alternative hits
# SA - chimeric read
samtools view -h /data/scratch/mpx469/stacks/ref-map/bwa/bwa-map-output/${INPUT_FILE}.mapped.sam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -h > samtools-output/${INPUT_FILE}.mapped.unique.sam

# test if all mapped reads are unique
if [ "`samtools view samtools-output/${INPUT_FILE}.mapped.unique.sam |  wc -l`" == "`samtools view samtools-output/${INPUT_FILE}.mapped.unique.sam | cut -f 1 | sort | uniq | wc -l`" ]; then
   echo "all good";
else
   echo "not good"
fi

# sort from name order into coordinate order
samtools sort -O bam -o samtools-output/${INPUT_FILE}.mapped.unique.sorted.bam -T ${INPUT_FILE}.tmp samtools-output/${INPUT_FILE}.mapped.unique.sam

# create index to allow viewing in igv
samtools index samtools-output/${INPUT_FILE}.mapped.unique.sorted.bam

