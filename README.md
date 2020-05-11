# Enset-tGBS
Enset tGBS methodology

This readme file details the methodolgy used in the analysis of ensete tGBS data. All work was completed on the QMUL apocrita cluster in the follwowing directory:
```
/data/scratch/mpx469/tGBS_enset_project
```

## Table of contents

[Import data from local hard drive into apocrita](#import-data-from-local-hard-drive-into-apocrita)

[Pre-processing of tGBS data](#pre-processing-of-tgbs-data)
   - [Create sample list to iterate through](#create-sample-list-to-iterate-through)
   - [Use trimmomatic to filter raw reads](#use-trimmomatic-to-filter-raw-reads)
   - [Use cutadapt to filter reads without restriction enzyme cutsites](#use-cutadapt-to-filter-reads-without-restriction-enzyme-cutsites)
   - [Get sample read counts for data2bio raw, trimmomatic and cutadapt](#get-sample-read-counts-for-data2bio-raw-trimmomatic-and-cutadapt)

[Stacks ref map pipeline](#stacks-ref-map-pipeline)
   - [BWA](#bwa)
   - [Samtools](#samtools)
   - [Gstacks](#gstacks)
   - [Populations](#populations)

[Post-processing of SNP data](#post-processing-of-snp-data)
   - [Filtering duplicated loci](#filtering-duplicated-loci)  
   - [Create blacklist for duplicated loci](#create-blacklist-for-duplicated-loci) 
   - [Repeat populations with blacklists](#repeat-populations-with-blacklists)


## Import data from local hard drive into apocrita
Set up directory for raw data

```
mkdir /data/scratch/mpx469/tGBS_enset_project/Data2Bio_final
```

Run in local terminal

```
rsync -avz --partial /drives/f/Genomic_data/Data2Bio_final/raw mpx469@login.hpc.qmul.ac.uk:/data/scratch/mpx469/tGBS_enset_project/Data2Bio_final
```
Set file permissions of Data2Bio directory

```
chmod -R u=rwx,g=r,o=r /data/scratch/mpx469/tGBS_enset_project/Data2Bio_final
```


## Pre-processing of tGBS data

### Create sample list to iterate through

```
# set dir
cd /data/scratch/mpx469/tGBS_enset_project
```

Import GBS_metadata in .csv format and create a sample list to iterate through

```
cut -f 2 -d "," tGBS_metadata.csv | tail -n +2 > sample_list.txt
```



### Use trimmomatic to filter raw reads


```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/trimmomatic
mkdir /data/scratch/mpx469/tGBS_enset_project/trimmomatic/trimmomatic_output

cd /data/scratch/mpx469/tGBS_enset_project/trimmomatic

qsub script_trimmomatic_array.sh

# all jobs should have run successfully
cat job_trimmomatic_array.o* | grep "done" -c
# should return 283
```



### Use cutadapt to filter reads without restriction enzyme cutsites

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/cutadapt
mkdir /data/scratch/mpx469/tGBS_enset_project/cutadapt/cutadapt_output

cd /data/scratch/mpx469/tGBS_enset_project/cutadapt

qsub script_cutadapt_array.sh

# all jobs should have run successfully
cat job_cutadapt_array.o* | grep done -c
# should return 283
```



### Use process radtags to truncate reads to uniform length

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/process_radtags

# will truncate reads at varying read lengths

# create output files
for i in `seq 70 10 120`; do  
   mkdir /data/scratch/mpx469/tGBS_enset_project/process_radtags/process_radtags_${i}_output
done

# submit jobs
for i in `seq 70 10 120`; do  
   qsub script_process_radtags_array_${i}.sh
done

# check output
for i in `seq 70 10 120`; do  
   cat job_process_radtags_${i}.o* | grep done -c
done
# should all equal 283
```


### map reads against the Bedadeti reference genome assembly

```
mkdir /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti

for i in `seq 70 10 120`; do  
   mkdir /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti/bwa_mem_${i}_output
   mkdir /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti/samtools_${i}_output
done

# index genome 
qsub script_bwa_index.sh

# map reads with bwa
for i in `seq 70 10 120`; do  
   qsub script_bwa_mem_${i}_array.sh
done

# check all jobs finished
for i in `seq 70 10 120`; do  
   cat job_bwa_mem_${i}_array.* | grep done -c 
done

# filter and sort with samtools
for i in `seq 70 10 120`; do  
   qsub script_samtools_${i}_array.sh
done

# check all jobs finished
for i in `seq 70 10 120`; do  
   cat job_samtools_${i}_array.* | grep done -c 
done

```



### assemble loci with gstacks 

```
mkdir /data/scratch/mpx469/tGBS_enset_project/gstacks

for i in `seq 70 10 120`; do  
   mkdir /data/scratch/mpx469/tGBS_enset_project/gstacks/gstacks_${i}_output
done

cd /data/scratch/mpx469/tGBS_enset_project/gstacks
```

Create popmap, filtering out samples classed as "Disease" or "NA" and submit script

```
Rscript write_popmap.R

qsub script_gstacks_array.sh
```


### call snps with populations

```
mkdir /data/scratch/mpx469/tGBS_enset_project/populations

for i in `seq 70 10 120`; do  
   mkdir /data/scratch/mpx469/tGBS_enset_project/populations/populations_${i}_single_snp_output
   mkdir /data/scratch/mpx469/tGBS_enset_project/populations/populations_${i}_all_snps_output
done

cd /data/scratch/mpx469/tGBS_enset_project/populations

qsub script_populations_all_snps_array.sh
qsub script_populations_single_snp_array.sh

grep "Populations is done" -c job_populations_*
```







### Get sample read counts for data2bio raw, trimmomatic and cutadapt

```
mkdir /data/scratch/mpx469/read-count
cd /data/scratch/mpx469/read-count

qsub script-read-number-count.sh
```

Get summary statisitics and create plots

```
module add R/3.6.1
Rscript Rscript-read-number-count-summary.R
```

**Table of summarys statistics** 
A full summary on an invididual basis can be found in the supplementary materials.

|         | raw (M) | trimmomatic (M) | removed by trimmomatic (%) | cutadapt (M)| removed by cutadapt (%) |
|---------|---------|-----------------|----------------------------|-------------|-------------------------|
| Total   |	832.11  | 639.73          | 0.23                       | 623.05      | 0.03                    |
| Mean	 | 3.20    | 2.46            | 0.23                       | 2.40        | 0.03                    |
| Min	    | 0.47    | 0.37            | 0.20                       | 0.36        | 0.01                    |
| Max	    | 11.80   | 9.18            | 0.32                       | 9.02        | 0.08                    |


**Read count histograms and boxplot**
![plot-read-count-histograms-boxplot](figures/plot-read-count-histograms-boxplot.png)


## Stacks ref map pipeline

### BWA

Create genome index and map reads to genome using BWA 

```
mkdir /data/scratch/mpx469/bwa/bwa-mem-output
mkdir /data/scratch/mpx469/bwa/bwa-mem-jobfiles
mkdir /data/scratch/mpx469/bwa/cutadapt-fq
cd /data/scratch/mpx469/bwa

# copy across genome
cp /data/SBCS-Ethiopia/databases/genomes/enset/GCA_000331365.3_Ensete_JungleSeeds_v3.0_genomic.fna.gz .

# gunzip
gunzip GCA_000331365.3_Ensete_JungleSeeds_v3.0_genomic.fna.gz

# index genome
qsub script-bwa-index.sh

# map trimmomatic output to reference
qsub script-bwa-mem.sh

# tidy up job files
mv job-bwa-mem.o* bwa-mem-jobfiles/

# should be 283
cat bwa-mem-jobfiles/job-bwa-mem.o* | grep -e "Real time" -c
```



### Samtools

Filter reads that are not mapped uniquely, flagged as alternative (XA:Z) or chimeric (SA:Z) alignments

```
mkdir /data/scratch/mpx469/samtools
mkdir /data/scratch/mpx469/samtoools/samtools-output
mkdir /data/scratch/mpx469/samtoools/samtools-job-files

cd /data/scratch/mpx469/samtools

qsub script-samtools.sh

# tidy up job files
mv job-samtools.o* samtools-job-files/

# should be 283
cat samtools-job-files/job-samtools.* | grep -e "all mapped reads are unique" -c

# get summary counts for mapped reads
qsub script-summary-flagstat-counts.sh

module add R/3.6.1
Rscript Rscript-summary-flagstat-counts.R
```

**Total and mapped reads from BWA and samtools**

![plot-flagstat](figures/plot-flagstat.png)


### Gstacks

```
mkdir /data/scratch/mpx469/STACKS
mkdir /data/scratch/mpx469/STACKS/gstacks
mkdir /data/scratch/mpx469/STACKS/gstacks/gstacks-together-output
mkdir /data/scratch/mpx469/STACKS/gstacks/gstacks-separate-output

cd /data/scratch/mpx469/STACKS/gstacks
```

Create popmap, filtering out samples classed as "Disease" or "NA". two pomap files are generate, one grouping all individuals as a single population, the other treating samples separately. 

```
Rscript Rscript-write-popmap-selection.R

 head popmap-selection-*
==> popmap-selection-separate.txt <==
EXS11ID000851.mapped.unique.sorted      1
EXS12ID000488.mapped.unique.sorted      2
EXS13ID000895.mapped.unique.sorted      3
EXS16ID000913.mapped.unique.sorted      4
EXS17ID000914.mapped.unique.sorted      5
EXS1ID000384.mapped.unique.sorted       6
EXS5ID000604.mapped.unique.sorted       7
EXS8ID000683.mapped.unique.sorted       8
P1EN003.mapped.unique.sorted    9
P1EN004.mapped.unique.sorted    10

==> popmap-selection-together.txt <==
EXS11ID000851.mapped.unique.sorted      1
EXS12ID000488.mapped.unique.sorted      1
EXS13ID000895.mapped.unique.sorted      1
EXS16ID000913.mapped.unique.sorted      1
EXS17ID000914.mapped.unique.sorted      1
EXS1ID000384.mapped.unique.sorted       1
EXS5ID000604.mapped.unique.sorted       1
EXS8ID000683.mapped.unique.sorted       1
P1EN003.mapped.unique.sorted    1
P1EN004.mapped.unique.sorted    1
```

Run gstacks

```
qsub script-gstacks-together.sh
qsub script-gstacks-separate.sh
```



### Populations

Run populatios using the popmaps created above

```
mkdir /data/scratch/mpx469/STACKS/populations

mkdir /data/scratch/mpx469/STACKS/populations/populations-separate-all-snps-R80-default-1
mkdir /data/scratch/mpx469/STACKS/populations/populations-separate-all-snps-R80-maf-het-1
mkdir /data/scratch/mpx469/STACKS/populations/populations-separate-single-snp-R80-default-1
mkdir /data/scratch/mpx469/STACKS/populations/populations-separate-single-snp-R80-maf-het-1
mkdir /data/scratch/mpx469/STACKS/populations/populations-together-all-snps-r80-default-1

cd /data/scratch/mpx469/STACKS/populations

qsub script-populations-separate-all-snps-R80-default-1.sh
qsub script-populations-separate-all-snps-R80-maf-het-1.sh
qsub script-populations-separate-single-snp-R80-default-1.sh
qsub script-populations-separate-single-snp-R80-maf-het-1.sh
qsub script-populations-together-all-snps-r80-default-1.sh

```


## Post-processing of SNP data

### Filtering duplicated loci 

Stacks assembles and defines loci with the same 5' startng point. A smll number of tGBS loci have different 5' start posiionts and appear in the vcf file as duplicate loci. An example of this is shown below, idenfited in IGV viewer.

** insert more recent image **

These loci need to be filtered out prior to analyses. 

#### Create blacklist for duplicated loci 

```
mkdir /data/scratch/mpx469/STACKS/blacklist-duplicates
mkdir /data/scratch/mpx469/STACKS/blacklist-duplicates/blacklist-separate-all-snps-R80-default
mkdir /data/scratch/mpx469/STACKS/blacklist-duplicates/blacklist-separate-all-snps-R80-maf-het
mkdir /data/scratch/mpx469/STACKS/blacklist-duplicates/blacklist-separate-single-snp-R80-default
mkdir /data/scratch/mpx469/STACKS/blacklist-duplicates/blacklist-separate-single-snp-R80-maf-het
mkdir /data/scratch/mpx469/STACKS/blacklist-duplicates/blacklist-together-all-snps-r80-default

bash script-blacklist-separate-all-snps-R80-default.sh
bash script-blacklist-separate-all-snps-R80-maf-het.sh
bash script-blacklist-separate-single-snp-R80-default.sh
bash script-blacklist-separate-single-snp-R80-maf-het.sh
bash script-blacklist-together-all-snps-r80-default.sh
```


#### Repeat populations with blacklists

```
cd /data/scratch/mpx469/STACKS/populations

mkdir /data/scratch/mpx469/STACKS/populations/populations-separate-all-snps-R80-default-2-blacklist
mkdir /data/scratch/mpx469/STACKS/populations/populations-separate-all-snps-R80-maf-het-2-blacklist
mkdir /data/scratch/mpx469/STACKS/populations/populations-separate-single-snp-R80-default-2-blacklist
mkdir /data/scratch/mpx469/STACKS/populations/populations-separate-single-snp-R80-maf-het-2-blacklist
mkdir /data/scratch/mpx469/STACKS/populations/populations-together-all-snps-r80-default-2-blacklist

qsub script-populations-separate-all-snps-R80-default-2-blacklist.sh
qsub script-populations-separate-all-snps-R80-maf-het-2-blacklist.sh
qsub script-populations-separate-single-snp-R80-default-2-blacklist.sh
qsub script-populations-separate-single-snp-R80-maf-het-2-blacklist.sh
qsub script-populations-together-all-snps-r80-default-2-blacklist.sh

```

Count snps in each dataset

```
grep ^# -v -c populations-*/populations.snps.vcf | sed -e 's/populations-//g' -e 's:/populations.snps.vcf::g'
separate-all-snps-R80-default-1:204436
separate-all-snps-R80-default-2-blacklist:204164
separate-all-snps-R80-maf-het-1:23707
separate-all-snps-R80-maf-het-2-blacklist:23698
separate-single-snp-R80-default-1:28729
separate-single-snp-R80-default-2-blacklist:28725
separate-single-snp-R80-maf-het-1:11193
separate-single-snp-R80-maf-het-2-blacklist:11189
together-all-snps-r80-default-1:204437
together-all-snps-r80-default-2-blacklist:204165
```

Copy populations output to shared space

```
cp -r /data/scratch/mpx469/STACKS/populations/ /data/SBCS-Ethiopia/tGBS_enset_project/
```







