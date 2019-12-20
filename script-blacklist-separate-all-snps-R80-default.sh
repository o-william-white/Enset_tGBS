
POPULATIONS_PATH="/data/scratch/mpx469/STACKS/populations/populations-separate-all-snps-R80-default-1/"
OUTPUT_PATH="blacklist-separate-all-snps-R80-default"

printf "\n Copying vcf and sumstats file to dir. \n"

# cp vcf file to dir
cp ${POPULATIONS_PATH}/populations.snps.vcf ${OUTPUT_PATH}/

# cp sumstats file to dir
cp ${POPULATIONS_PATH}/populations.sumstats.tsv ${OUTPUT_PATH}/

printf "\n Counting total, duplicate and unique snps in vcf file \n"

# total number of snps in vcf file 
SNPS_TOTAL=`grep -e ^# -v ${OUTPUT_PATH}/populations.snps.vcf | wc -l`

# number of snps with the same chromosome and position
# -D prints all duplicate lines
SNPS_DUPS=`grep -e ^# -v ${OUTPUT_PATH}/populations.snps.vcf | awk ' { print $1"\t"$2 } ' | sort | uniq -D | wc -l`

# number of snps after filtering those with identical chromosome and position
# -u only print unique lines
SNPS_UNIQ=`grep -e ^# -v ${OUTPUT_PATH}/populations.snps.vcf | awk ' { print $1"\t"$2 } ' | sort | uniq -u | wc -l`

printf "\n Writing file with duplicated snps called out-duplicated-sites.txt \n"
printf "   Col 1 = number of occurences \n"
printf "   Col 2 = chr. \n"
printf "   Col 3 = bp \n"

# create file with duplicated snps and count
# uniq -d only print duplicate lines, one for each group  
#      -c prefix lines by the number of occurrences
grep -e ^# -v ${OUTPUT_PATH}/populations.snps.vcf | awk ' { print $1"\t"$2 } ' | sort | uniq -dc > ${OUTPUT_PATH}/out-duplicated-sites.txt


# need to link chromosome name and bp with locus id for blacklist

# get locus id, chromosome name and bp from sumstats file
printf "\n Writing file called out-sumstats-info.txt \n"
printf "   Col 1 = locus id \n"
printf "   Col 2 = chr. \n"
printf "   Col 3 = bp \n"

grep -e "^#" -v ${OUTPUT_PATH}/populations.sumstats.tsv | awk ' { print $1"\t"$2"\t"$3 } ' | uniq > ${OUTPUT_PATH}/out-sumstats-info.txt

# join snp chromosome name and bp to locus number in R
printf "\n Joining chromosome name and bp to locus number in R \n"
printf " Writing output to a file called out-duplicated-sites-and-sumstats-info.txt \n"
printf "   Col 1 = number of occurences \n"
printf "   Col 2 = chr. \n"
printf "   Col 3 = bp \n"
printf "   Col 4 = locus \n\n"

module add R/3.6.1

# run R script
# note path must be specified
Rscript   Rscript-join.R   ${OUTPUT_PATH}

# write blacklist
printf "\n Writing blacklist \n"
cut -f 4 ${OUTPUT_PATH}/out-duplicated-sites-and-sumstats-info.txt | sort | uniq > ${OUTPUT_PATH}/blacklist.txt

LOCI=`wc -l ${OUTPUT_PATH}/blacklist.txt | cut -f 1 -d " "` 

printf "\n Summary \n"

printf "\n No. of snps in vcf file:                                  ${SNPS_TOTAL}"
printf "\n No. with the same chr. and pos. :                         ${SNPS_DUPS}"
printf "\n No. after filtering those with identical chr. and pos. :  ${SNPS_UNIQ}"
printf "\n No. of loci added to blacklist                            ${LOCI} \n  "  

printf "\n Summary written to ${OUTPUT_PATH}/logfile \n"

echo "No. of snps in vcf file:                                  ${SNPS_TOTAL}" > ${OUTPUT_PATH}/logfile
echo "No. with the same chr. and pos. :                         ${SNPS_DUPS}" >> ${OUTPUT_PATH}/logfile
echo "No. after filtering those with identical chr. and pos. :  ${SNPS_UNIQ}" >> ${OUTPUT_PATH}/logfile
echo "No. of loci added to blacklist                            ${LOCI}     " >> ${OUTPUT_PATH}/logfile

printf "\n Done \n"

