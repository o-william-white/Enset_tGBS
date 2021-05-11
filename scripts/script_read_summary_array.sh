#!/bin/sh
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -t 1-283
#$ -tc 36
#$ -N job_read_summary_array

SAMPLE=$(sed -n "${SGE_TASK_ID}p" /data/scratch/mpx469/tGBS_enset_project/sample_list.txt)

module load samtools 

# set up output prior to running script
# for LENGTH in `seq 70 10 120`; do
#   echo SAMPLE,RAW,TRIMMOMATIC,CUTADAPT,CUTADAPT_BP,PROCESS_RADTAGS,PROCESS_RADTAGS_BP,BWA_UNMAPPED,BWA_MAPPED,BWA_XA_SA,SAM_UNMAPPED,SAM_MAPPED > summary_$LENGTH.csv
# done

RAW=$(zcat /data/scratch/mpx469/tGBS_enset_project/Data2Bio_final/raw/${SAMPLE}.digested.fq.gz | echo $((`wc -l`/4)))
TRIMMOMATIC=$(zcat /data/scratch/mpx469/tGBS_enset_project/trimmomatic/trimmomatic_output/${SAMPLE}.trimmomatic.fq.gz | echo $((`wc -l`/4)))
CUTADAPT=$(zcat /data/scratch/mpx469/tGBS_enset_project/cutadapt/cutadapt_output/${SAMPLE}.cutadapt.with.cutsites.fq.gz | echo $((`wc -l`/4)))
CUTADAPT_BP=$(zcat /data/scratch/mpx469/tGBS_enset_project/cutadapt/cutadapt_output/${SAMPLE}.cutadapt.with.cutsites.fq.gz | paste - - - - | cut -f 2 | tr -d '\n' | wc -c)

for LENGTH in `seq 70 10 120`; do
   PROCESS_RADTAGS=$(zcat /data/scratch/mpx469/tGBS_enset_project/process_radtags/process_radtags_${LENGTH}_output/${SAMPLE}.cutadapt.with.cutsites.fq.gz | echo $((`wc -l`/4)))
   PROCESS_RADTAGS_BP=$(zcat /data/scratch/mpx469/tGBS_enset_project/process_radtags/process_radtags_${LENGTH}_output/${SAMPLE}.cutadapt.with.cutsites.fq.gz | paste - - - - | cut -f 2 | tr -d '\n' | wc -c)
   BWA_UNMAPPED=$(samtools view -f 0x4 /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti/bwa_mem_${LENGTH}_output/${SAMPLE}.bam  | wc -l)
   BWA_MAPPED=$(samtools view -F 0x4 /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti/bwa_mem_${LENGTH}_output/${SAMPLE}.bam  | wc -l)
   BWA_XA_SA=$(samtools view /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti/bwa_mem_${LENGTH}_output/${SAMPLE}.bam | grep -e 'XA:Z:' -e 'SA:Z:' -c)
   SAM_UNMAPPED=$(samtools view -f 0x4 /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti/samtools_${LENGTH}_output/${SAMPLE}.bam  | wc -l)
   SAM_MAPPED=$(samtools view -F 0x4 /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti/samtools_${LENGTH}_output/${SAMPLE}.bam  | wc -l)
   echo $SAMPLE,$RAW,$TRIMMOMATIC,$CUTADAPT,$CUTADAPT_BP,$PROCESS_RADTAGS,$PROCESS_RADTAGS_BP,$BWA_UNMAPPED,$BWA_MAPPED,$BWA_XA_SA,$SAM_UNMAPPED,$SAM_MAPPED >> summary_$LENGTH.csv
done
   
echo $SAMPLE complete

