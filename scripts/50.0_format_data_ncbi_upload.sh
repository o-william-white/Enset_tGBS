#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=6:0:0
#$ -l h_vmem=1G
#$ -j y
#$ -N job_format_data_ncbi_upload

mkdir -p Data2Bio_final/ncbi_upload

ls -1 Data2Bio_final/raw/*.fq.gz | while read LINE; do 
   
   BASENAME=$(basename $LINE | sed -e 's/.fq.gz/.fastq.gz/g' -e 's/.digested/_digested/g')
   
   echo copying $LINE
   cp $LINE Data2Bio_final/ncbi_upload/${BASENAME}
   
   echo checking md5sum
   md5sum $LINE > Data2Bio_final/ncbi_upload/${BASENAME}.md5
   md5sum -c Data2Bio_final/ncbi_upload/${BASENAME}.md5

done

echo Complete

