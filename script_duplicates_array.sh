#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N job_duplicates_array
#$ -t 1-12

LENGTH=$(sed -n "${SGE_TASK_ID}p" input_args | cut -f 1 -d " ")
DATASET=$(sed -n "${SGE_TASK_ID}p" input_args | cut -f 2 -d " ")

# get populations and output paths
POPULATIONS_PATH="/data/scratch/mpx469/tGBS_enset_project/populations/populations_${LENGTH}_${DATASET}_output/"
OUTPUT_PATH="duplicates_${LENGTH}_${DATASET}_output"

# mk output dir
mkdir -p ${OUTPUT_PATH}

# cp vcf file to dir
cp ${POPULATIONS_PATH}/populations.snps.vcf ${OUTPUT_PATH}/

# cp sumstats file to dir
cp ${POPULATIONS_PATH}/populations.sumstats.tsv ${OUTPUT_PATH}/

# total number of snps in vcf file
SNPS_TOTAL=`grep -e ^# -v ${OUTPUT_PATH}/populations.snps.vcf | wc -l`

# number of snps with the same chromosome and position
# -D prints all duplicate lines
SNPS_DUPS=`grep -e ^# -v ${OUTPUT_PATH}/populations.snps.vcf | awk ' { print $1"\t"$2 } ' | sort | uniq -D | wc -l`

# number of snps after filtering those with identical chromosome and position
# -u only print unique lines
SNPS_UNIQ=`grep -e ^# -v ${OUTPUT_PATH}/populations.snps.vcf | awk ' { print $1"\t"$2 } ' | sort | uniq -u | wc -l`


# create file with duplicated snps and count
# uniq -d only print duplicate lines, one for each group
#      -c prefix lines by the number of occurrences
# col 1 = number of occurences
# col 2 = chr.
# col 3 = bp
grep -e ^# -v ${OUTPUT_PATH}/populations.snps.vcf | awk ' { print $1"\t"$2 } ' | sort | uniq -dc > ${OUTPUT_PATH}/duplicates


# need to link chromosome name and bp with locus id for blacklist

# get locus id, chromosome name and bp from sumstats file
# col 1 = locus id
# col 2 = chr.
# col 3 = bp
grep -e "^#" -v ${OUTPUT_PATH}/populations.sumstats.tsv | awk ' { print $1"\t"$2"\t"$3 } ' | uniq > ${OUTPUT_PATH}/loci_info

# join snp chromosome name and bp to locus number in R
# output file called duplicates_loci'
# col 1 = number of occurences
# col 2 = chr.
# col 3 = bp
# col 4 = locus

module add R/3.6.1

# run R script
# note path must be specified
Rscript Rscript_join_loci_numbers.R ${OUTPUT_PATH}

# write blacklist
cut -f 4 ${OUTPUT_PATH}/duplicates_loci | sort | uniq > blacklist_duplicates_${LENGTH}_${DATASET}

# number of loci in blacklist
BLACKLIST=`wc -l blacklist_duplicates_${LENGTH}_${DATASET} | cut -f 1 -d " "`

echo "total      ${SNPS_TOTAL}" > ${OUTPUT_PATH}/summary
echo "duplicate  ${SNPS_DUPS}" >> ${OUTPUT_PATH}/summary
echo "remaining  ${SNPS_UNIQ}" >> ${OUTPUT_PATH}/summary
echo "blacklist  ${BLACKLIST}" >> ${OUTPUT_PATH}/summary

printf "\n summary written to ${OUTPUT_PATH}/summary \n"

# remove unnecessary files
rm ${OUTPUT_PATH}/populations.snps.vcf
rm ${OUTPUT_PATH}/populations.sumstats.tsv
