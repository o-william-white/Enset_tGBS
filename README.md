# ensete-tGBS
Ensete tGBS methodology

This readme file details the methodolgy used in the analysis of ensete tGBS data. All work was completed in QMUL apocrita cluster.

### Import data from local hard drive into apocrita
Set up directory for raw data

`mkdir /data/scratch/mpx469/Data2Bio_final`

Run in local terminal

```
rsync -avz --partial /drives/f/Genomic_data/Data2Bio_final/raw mpx469@login.hpc.qmul.ac.uk:/data/scratch/mpx469/Data2Bio_final
rsync -avz --partial /drives/f/Genomic_data/Data2Bio_final/trimmed mpx469@login.hpc.qmul.ac.uk:/data/scratch/mpx469/Data2Bio_final
rsync -avz --partial /drives/f/Genomic_data/Data2Bio_final/genome mpx469@login.hpc.qmul.ac.uk:/data/scratch/mpx469/Data2Bio_final
rsync -avz --partial /drives/f/Genomic_data/Data2Bio_final/alignment.BAM mpx469@login.hpc.qmul.ac.uk:/data/scratch/mpx469/Data2Bio_final
```
Set file permissions of Data2Bio directory

```
chmod -R u=rwx,g=r,o=r Data2Bio_final
```

Import sample meta data into /data/scratch/mpx469 in .csv and .txt formats

```
head GBS_metadata.txt
SEQUENCE_ID     SAMPLE_ID       POPULATION      TYPE    sample_id       Landrace_2      _latitude       _longitude      Unique_code     Sample_date
EXS11ID000851.digested.trimmed.fq.gz    EXS11ID000851   1       Domestic        851     Chamo   6.191268003     37.57450757     031c0c90-0ac4-4f68-8f69-d960eaecb705    2018-10-18
EXS12ID000488.digested.trimmed.fq.gz    EXS12ID000488   1       Domestic        488     Worsaife        6.049165501     37.22941625     5a728129-b1e3-4781-9890-756b87133039    2018-10-20
EXS13ID000895.digested.trimmed.fq.gz    EXS13ID000895   1       Domestic        895     Arooko  7.933571263     36.51323982     c6b997b2-95b1-4b8b-89e3-a15fbeb4303e    2018-10-24
EXS16ID000913.digested.trimmed.fq.gz    EXS16ID000913   1       Domestic        913     Unknown_green   9.037241562     37.42883203     d873d0e5-70ef-4e53-a47f-144583374165    2018-10-25
EXS17ID000914.digested.trimmed.fq.gz    EXS17ID000914   1       Domestic        914     Unknown_red     9.037225949     37.4289773      d873d0e5-70ef-4e53-a47f-144583374165    2018-10-25
EXS1ID000384.digested.trimmed.fq.gz     EXS1ID000384    1       Domestic        384     Ganticha        6.774082442     38.43834313     835c3a82-4e2d-4851-bc81-edafefc6b651    2018-04-02
EXS5ID000604.digested.trimmed.fq.gz     EXS5ID000604    1       Domestic        604     Toracho 6.137371431     38.1996072      8b3cb20c-10fd-41c4-8cf4-0316c66bde65    2018-04-04
EXS8ID000683.digested.trimmed.fq.gz     EXS8ID000683    1       Domestic        683     Wanade  6.778954922     37.76874955     7596b9e7-27ca-4a7c-987f-426956b161e2    2018-04-10
P1EN003.digested.trimmed.fq.gz  P1EN003 1       Domestic        162     Achachet        8.491123519     38.01776645     0f21ebee-c784-4174-9e60-a046b3f5ab97    2018-01-27
```



### Create sample list to iterate through

```
# set dir
cd /data/scratch/mpx469

# get list from metadata
cut -f 2 GBS_metadata.txt | tail -n +2 > sample-list.txt

```



### Use trimmomatic to filter raw reads

```
# set dir
mkdir /data/scratch/mpx469/trimmomatic
mkdir /data/scratch/mpx469/trimmomatic/trimmomatic-output
cd /data/scratch/mpx469/trimmomatic

qsub script-trimmomatic-array.sh
```

Note LEADING paramter not used, in attempt to preserve the 5' end of the GBS loci which are used to align stacks

```
# tidy up jobfiles
mkdir trimmomatic-output/jobfiles
cp job-trimmomatic-array.o* trimmomatic-output/jobfiles/

# all jobs should have run successfully
cat trimmomatic-output/jobfiles/job-trimmomatic-array.o* | grep -e "TrimmomaticSE: Completed successfully" -c
# should return 283
```



### Get sample read counts for data2bio raw, data3bio trimmed and trimmomatic

```
mkdir /data/scratch/mpx469/read-count
cd /data/scratch/mpx469/read-count

qsub script-read-number-count.sh
```

Create plots

```
module add R
Rscript Rscript-plot-read-number-count.R
```

![GitHub Logo](/figure/plot-read-count-histograms.png)
Format: ![Alt Text](url)


### Read length distributions

Ion proton sequencing has variable read lengths. Quantify read lengths in a given sample, raw and trimmed. Code adapted from https://www.biostars.org/p/72433/#72441

```
mkdir /data/scratch/mpx469/read-length-distribution
cd /data/scratch/mpx469/read-length-distribution

# compare read length distributions of a single sample for the different trimming options
mkdir /data/scratch/mpx469/read-length-distribution/single-sample-comparison
cd /data/scratch/mpx469/read-length-distribution/single-sample-comparison/

# get read lengths
zcat /data/scratch/mpx469/Data2Bio_final/raw/EXS11ID000851.digested.fq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > read_length.EXS11ID000851.digested.txt
zcat /data/scratch/mpx469/Data2Bio_final/trimmed/EXS11ID000851.digested.trimmed.fq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > read_length.EXS11ID000851.digested.trimmed.txt
zcat /data/scratch/mpx469/trimmomatic/trimmomatic-output/EXS11ID000851.digested.trimmomatic.fq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > read_length.EXS11ID000851.digested.trimmomatic.txt

# column 1 is the number of sequnces
# column 2 is the length category
```

Create plot

```
module add R

cat Rscript-plot-single-sample-comparison.R
# read in input data
reads.raw         <- read.table(file="read_length.EXS11ID000851.digested.txt", header=FALSE)
reads.trim        <- read.table(file="read_length.EXS11ID000851.digested.trimmed.txt", header=FALSE)
reads.trimmomatic <- read.table(file="read_length.EXS11ID000851.digested.trimmomatic.txt", header=FALSE)

# create line plot
pdf(file="plot.read_length.EXS11ID000851.digested.pdf")
plot(reads.raw$V2,reads.raw$V1,type="l",xlab="read length",ylab="occurences",col="blue", main="EXS11ID000851")
lines(reads.trim$V2,reads.trim$V1, col="red")
lines(reads.trimmomatic$V2,reads.trimmomatic$V1, col="green")
legend(min(reads.raw$V2), max(reads.raw$V1), legend=c("raw", "trimmed", "trimmomatic"), col=c("blue", "red", "green"), lty=c(1,1,1), ncol=1)
dev.off()

q(save="no")

Rscript Rscript-plot-single-sample-comparison.R
```

Get read length distributions for all trimmomatic samples

```
mkdir /data/scratch/mpx469/read-length-distribution/read-length-distribution-trimmomatic
mkdir /data/scratch/mpx469/read-length-distribution/read-length-distribution-trimmomatic/output-read-length-distribution-trimmomatic

cd /data/scratch/mpx469/read-length-distribution/read-length-distribution-trimmomatic/

cat script-read-length-distribution-trimmomatic.sh
#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=2:0:0
#$ -cwd
#$ -j y
#$ -N job-read-length-distribution-trimmomatic

cat /data/scratch/mpx469/sample-list.txt | while read i; do
   zcat /data/scratch/mpx469/trimmomatic/trimmomatic-output/$i'.digested.trimmomatic.fq.gz' | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > output-read-length-distribution-trimmomatic/read_length.$i'.digested.trimmomatic.txt'
done

qsub script-read-length-distribution-trimmomatic.sh
```
Create plot 

```
module add R/3.6.1
R

# install.packages('dplyr',   lib="/data/home/mpx469/software/R/3.6.1/", repos = 'https://cloud.r-project.org')
# install.packages('ggplot2', lib="/data/home/mpx469/software/R/3.6.1/", repos = 'https://cloud.r-project.org')

# set lib path
.libPaths("/data/home/mpx469/software/R/3.6.1/")

# import libraries
library(dplyr)
library(ggplot2)

input.files <- list.files(path="output-read-length-distribution-trimmomatic/", pattern="read_length*", full.names=TRUE)

input.list <- lapply(input.files, read.table)

head(input.list[[1]])
head(input.list[[2]])

# V1 = number of reads
# V2 = read length

# plot read length distribution of all input files
pdf("plot-read-dristribution-trimmomatic.pdf")
plot(0,0,xlim = c(0,250),ylim = c(0, 200000), type = "n", xlab="Read length", ylab="Number of reads")
lapply(input.list, function(x) lines(x$V2, x$V1, col = "blue"))
dev.off()

my_func <- function(x) {
  x %>%
    summarise(cut60  = sum(x[x$V2>=60 ,"V1"]),
              cut80  = sum(x[x$V2>=80 ,"V1"]),
              cut100 = sum(x[x$V2>=100,"V1"]),
              cut120 = sum(x[x$V2>=120,"V1"]),
              cut140 = sum(x[x$V2>=140,"V1"]),
              cut160 = sum(x[x$V2>=160,"V1"])) 
} 

read.length.number <- lapply(input.list, my_func)

# plot the number of reads after truncating read lengths
pdf("plot-read-numbers-after-truncating-lengths.pdf")
plot(0,0,xlim = c(60,160),ylim = c(0, 10000000), type = "n", xlab="Read length", ylab="Number of reads")
for (i in 1:length(read.length.number)) {
  lines(seq(60,160,by=20), read.length.number[[i]], col = "blue")
}
dev.off()

q(save="no")

# overcrowded plot
# is there a way to average across?
# perhaps more importantly, should we drop certain samples that have minimum number of samples?
```
 


# stacks pipeline reference guided - BWA

Create genome index and map reads to genome using BWA 

```
mkdir /data/scratch/mpx469/stacks/ref-map/
mkdir /data/scratch/mpx469/stacks/ref-map/bwa
cd /data/scratch/mpx469/stacks/ref-map/bwa

# copy across genome
cp /data/SBCS-Ethiopia/databases/genomes/enset/GCA_000331365.3_Ensete_JungleSeeds_v3.0_genomic.fna.gz .

gunzip GCA_000331365.3_Ensete_JungleSeeds_v3.0_genomic.fna.gz

# index genome
cat script-bwa-index.sh
#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G
#$ -N job-bwa-index

module add bwa

# index the reference genome
bwa index GCA_000331365.3_Ensete_JungleSeeds_v3.0_genomic.fna

qsub script-bwa-index.sh
```

Map trimmomatic output to reference

```
mkdir bwa-map-output
mkdir trimmomatic-fq

cat script-bwa-mem.sh
#!/bin/sh
#$ -cwd
#$ -j y
#$ -l h_rt=2:0:0
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -t 1-283
#$ -tc 283
#$ -N job-bwa-map

module add bwa

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" /data/scratch/mpx469/sample-list.txt)

cp /data/scratch/mpx469/trimmomatic/trimmomatic-output/${INPUT_FILE}.digested.trimmomatic.fq.gz ./trimmomatic-fq/

gunzip ./trimmomatic-fq/${INPUT_FILE}.digested.trimmomatic.fq.gz

bwa mem GCA_000331365.3_Ensete_JungleSeeds_v3.0_genomic.fna ./trimmomatic-fq/${INPUT_FILE}.digested.trimmomatic.fq > ${INPUT_FILE}.mapped.sam

qsub script-bwa-mem.sh

# tidy up job files
mkdir bwa-map-output/job-files
mv job-bwa-map.o* bwa-map-output/job-files/

# should be 283
cat bwa-map-output/job-files/job-bwa-map.o* | grep -e "Real time" -c
```



### stacks pipeline reference guided - samtools

```
mkdir /data/scratch/mpx469/stacks/ref-map/samtools
cd /data/scratch/mpx469/stacks/ref-map/samtools

mkdir samtools-output

cat script-samtools-filter-sort-index.sh
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

qsub script-samtools-filter-sort-index.sh


mkdir samtools-output/job-files
mv job-samtools.o* samtools-output/job-files/

cat samtools-output/job-files/job-samtools.* | grep -e "all good" -c
283
```

How many reads filtered and mapped 

```
cat script-samtools-flagstat.sh
#!/bin/sh
#$ -cwd
#$ -j y
#$ -l h_rt=1:0:0
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -N job-samtools-flagstat
#$ -t 1-283
#$ -tc 283

module add samtools

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" /data/scratch/mpx469/sample-list.txt)

samtools flagstat /data/scratch/mpx469/stacks/ref-map/bwa/bwa-map-output/${INPUT_FILE}.mapped.sam > samtools-output/${INPUT_FILE}.flagstat.txt

samtools flagstat /data/scratch/mpx469/stacks/ref-map/samtools/samtools-output/${INPUT_FILE}.mapped.unique.sam > samtools-output/${INPUT_FILE}.flagstat.filtered.txt

qsub script-samtools-flagstat.sh

# job files not needed or useful
rm job-samtools-flagstat.o*

mkdir flagstat-plot

cat /data/scratch/mpx469/sample-list.txt | while read i; do  
   head -n 1 samtools-output/${i}.flagstat.txt | cut -f 1 -d " "; 
done > flagstat-plot/reads-total.txt

cat /data/scratch/mpx469/sample-list.txt | while read i; do  
   head -n 1 samtools-output/${i}.flagstat.filtered.txt | cut -f 1 -d " "; 
done > flagstat-plot/reads-unique.txt

cat /data/scratch/mpx469/sample-list.txt | while read i; do  
   awk 'NR==5' samtools-output/${i}.flagstat.filtered.txt | cut -f 1 -d " "; 
done > flagstat-plot/reads-unique-mapped.txt

module add R
R

# read in data
total <- read.table("flagstat-plot/reads-total.txt", header=FALSE)
unique <- read.table("flagstat-plot/reads-unique.txt", header=FALSE)
unique.mapped <- read.table("flagstat-plot/reads-unique-mapped.txt", header=FALSE)

# rbind data in long format
df <- rbind(total, unique, unique.mapped)

# add class info
df <- cbind(df, c(rep("total", nrow(total)), rep("unique", nrow(unique)), rep("unique mapped", nrow(unique))))

colnames(df) <- c("count", "class")

pdf(file="flagstat-plot/plot-flagstat.pdf")
boxplot(count ~ class, df, xlab = "", sub=paste("Mean proportion of uniquely mapped reads = ", signif(mean(unique.mapped$V1 / unique$V1), digits=2), "%", sep =""))
dev.off()

q(save="no")
```



### stacks pipeline reference guided - gstacks

```
mkdir /data/scratch/mpx469/stacks/ref-map/gstacks
cd /data/scratch/mpx469/stacks/ref-map/gstacks

# run gstacks after removing unneccessary samples (Disease or NA)
# run as all sample lumped together and all samples treated separately


module add R
R

sample.metadata <- read.csv("/data/scratch/mpx469/GBS_metadata.csv", header=TRUE, na.strings=NULL)

# filter out disease and NA samples
sample.metadata <- sample.metadata[-which(sample.metadata$TYPE == "Disease" | sample.metadata$TYPE == "NA"),]

samples <- paste(sample.metadata$SAMPLE_ID, ".mapped.unique.sorted", sep="")

population.together <- rep(1, length(samples))
population.separate <- 1:length(samples)

popmap.together <- cbind(samples, population.together)
popmap.separate <- cbind(samples, population.separate)

write.table(file="popmap-selection-together.txt", popmap.together, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(file="popmap-selection-separate.txt", popmap.separate, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

q(save="no")

mkdir gstacks-together-output
mkdir gstacks-separate-output


cat script-gstacks-together.sh
#!/bin/sh
#$ -cwd
#$ -j y
#$ -pe smp 12
#$ -l h_rt=12:0:0
#$ -l h_vmem=2G
#$ -N job-gstacks-together

module load use.dev
module add stacks/2.41

gstacks \
   -I /data/scratch/mpx469/stacks/ref-map/samtools/samtools-output/ \
   -M popmap-selection-together.txt \
   -O gstacks-together-output \
   -t 12

cat script-gstacks-separate.sh
#!/bin/sh
#$ -cwd
#$ -j y
#$ -pe smp 12
#$ -l h_rt=12:0:0
#$ -l h_vmem=2G
#$ -N job-gstacks-separate

module load use.dev
module add stacks/2.41

gstacks \
   -I /data/scratch/mpx469/stacks/ref-map/samtools/samtools-output/ \
   -M popmap-selection-separate.txt \
   -O gstacks-separate-output \
   -t 12

qsub script-gstacks-together.sh
qsub script-gstacks-separate.sh
```



### stacks pipeline reference guided - populations

Run populations after removing unneccessary samples (Disease or NA)

Also run and as all sample lumped together as a single populaytion and and all samples treated separately

```
mkdir /data/scratch/mpx469/stacks/ref-map/populations
cd /data/scratch/mpx469/stacks/ref-map/populations

mkdir populations-together-output
mkdir populations-separate-output

cat script-populations-together.sh
#!/bin/bash
#$ -pe smp 12
#$ -l h_vmem=4G
#$ -l h_rt=2:0:0
#$ -cwd
#$ -j y
#$ -N job-populations-together

module load use.dev
module add stacks/2.41

populations \
   -P /data/scratch/mpx469/stacks/ref-map/gstacks/gstacks-together-output/ \
   -O populations-together-output \
   -M /data/scratch/mpx469/stacks/ref-map/gstacks/popmap-selection-together.txt \
   -r 0.4 \
   --vcf --plink --structure \
   -t 12

cat script-populations-separate.sh
#!/bin/bash
#$ -pe smp 12
#$ -l h_vmem=24G
#$ -l h_rt=10:0:0
#$ -l node_type=sm
#$ -l highmem
#$ -cwd
#$ -j y
#$ -N job-populations-separate

module load use.dev
module add stacks/2.41

populations \
   -P /data/scratch/mpx469/stacks/ref-map/gstacks/gstacks-separate-output/ \
   -O populations-separate-output \
   -M /data/scratch/mpx469/stacks/ref-map/gstacks/popmap-selection-separate.txt \
   -R 0.4 \
   --write-single-snp \
   --phylip --phylip-var --phylip-var-all --vcf \
   -t 12
```



### identify duplicated loci 

```
# how many mapped loci in total
grep -e ^# -v populations.snps.vcf | wc -l
576099

# how many mapped loci after filtering those with identical chromosome and position
# nout -u call needed to only print unique lines
grep -e ^# -v populations.snps.vcf | awk ' { print $1"\t"$2 } ' | sort | uniq -u | wc -l
570318

# how many lines represent duplicated loci
expr 576099 - 570318
5781

# check the number of duplicated lines
grep -e ^# -v populations.snps.vcf | awk ' { print $1"\t"$2 } ' | sort | uniq -D | wc -l
5781


# grep out info lines begining with #
# select first and second lines with chomosome and position
# sort
# uniq - filter for duplicate lines
#      -d only print duplicate lines, one for each group  
#      -c prefix lines by the number of occurrences
grep -e ^# -v populations.snps.vcf | awk ' { print $1"\t"$2 } ' | sort | uniq -dc > duplicated-sites.txt


module add R
R

x <- read.table("duplicated-sites.txt", header=FALSE)

head(x)
#  V1             V2    V3
#1  2 AMZH03000010.1 31526
#2  2 AMZH03000010.1 31541
#3  2 AMZH03000010.1 31546
#4  2 AMZH03000010.1 31548
#5  2 AMZH03000010.1 31555
#6  2 AMZH03000010.1 31564

sum(x$V1)
#[1] 5781

# add an empty column 
# I would like the to be a tab on the end of each line
x$V4 <- rep("", nrow(x))

write.table(x[,-1], "duplicated-sites-edit.txt", sep ="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

q(save="no")

grep -f duplicated-sites-edit.txt -c populations.snps.vcf
#5781

grep -f duplicated-sites-edit.txt -v populations.snps.vcf > populations.snps.duplicates.removed.vcf


grep -e ^# -v populations.snps.duplicates.removed.vcf | wc -l
#570318




wc -l duplicated-sites.txt
2871 duplicated-sites.txt

# what is the maximum number of duplicates
cut -f 7 -d " " duplicated-sites.txt | sort | uniq
2
3

grep -e ^# -v populations.snps.vcf | awk ' { print $1"\t"$2 } ' | sort | uniq -D > duplicated-sites-all.txt
```





