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
   - [Use process radtags to truncate reads to uniform length](#use-process-radtags-to-truncate-reads-to-uniform-length)
   - [map reads against the Bedadeti reference genome assembly](#map-reads-against-the-Bedadeti-reference-genome-assembly)
   - [assemble loci with gstacks](#assemble-loci-with-gstacks)
   - [call snps with populations](#call-snps-with-populations)
   - [blacklist contaminants and paralogues](#blacklist-contaminants-and-paralogues)
   - [](#)
   - [](#)


## Import data from local hard drive into apocrita

Set up directory for raw data
```
mkdir /data/scratch/mpx469/tGBS_enset_project/Data2Bio_final
```

Run in local terminal
```
rsync -avz --partial /drives/f/Genomic_data/Data2Bio_final/raw mpx469@login.hpc.qmul.ac.uk:/data/scratch/mpx469/tGBS_enset_project/Data2Bio_final
```

Set file permissions of Data2Bio directory in apocrita
```
chmod -R u=rwx,g=r,o=r /data/scratch/mpx469/tGBS_enset_project/Data2Bio_final
```

<br/>
<div align="right">
    <b><a href="#----">↥ back to top</a></b>
</div>
<br/>

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
```



### Use cutadapt to filter reads without restriction enzyme cutsites

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/cutadapt
mkdir /data/scratch/mpx469/tGBS_enset_project/cutadapt/cutadapt_output

cd /data/scratch/mpx469/tGBS_enset_project/cutadapt

qsub script_cutadapt_array.sh
```



### Use process radtags to truncate reads to uniform length

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/process_radtags

# will truncate reads at varying read lengths
# 60 bp to 120 bp in steps of 10

# create output files
for l in `seq 70 10 120`; do  
   mkdir /data/scratch/mpx469/tGBS_enset_project/process_radtags/process_radtags_${l}_output
done

# create input args for array script
for l in `seq 70 10 120`; do 
   cat /data/scratch/mpx469/tGBS_enset_project/sample_list.txt | while read i; do
      echo $l $i >> input_args
   done	  
done

qsub script_process_radtags_array.sh
```



### map reads against the Bedadeti reference genome assembly

```
mkdir /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti

for l in `seq 70 10 120`; do  
   mkdir /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti/bwa_mem_${l}_output
   mkdir /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti/samtools_${l}_output
done

cd /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti

# index genome 
qsub script_bwa_index.sh

# cp input_args to dir
cp /data/scratch/mpx469/tGBS_enset_project/process_radtags/input_args .

qsub script_map_reads_bedadeti_array.sh

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

for l in `seq 70 10 120`; do  
   mkdir /data/scratch/mpx469/tGBS_enset_project/populations/populations_${l}_single_snp_output
   mkdir /data/scratch/mpx469/tGBS_enset_project/populations/populations_${l}_all_snps_output
done

cd /data/scratch/mpx469/tGBS_enset_project/populations

qsub script_populations_all_snps_array.sh
qsub script_populations_single_snp_array.sh

# get summary stats for each assembly
echo -e 'loci sites filtered variant' > summary_all_snps
echo -e 'loci sites filtered variant' > summary_single_snp

for l in `seq 70 10 120`; do 
    echo ${l} $(grep Kept populations_${l}_all_snps_output/populations.log | cut -f 2,6,8,14 -d " ") >> summary_all_snps
    echo ${l} $(grep Kept populations_${l}_single_snp_output/populations.log | cut -f 2,6,8,14 -d " ") >> summary_single_snp
done

Rscript plot_summary.R

```

<br/>
<div align="right">
    <b><a href="#----">↥ back to top</a></b>
</div>
<br/>

## blacklist contaminants and paralogues

```
mkdir /data/scratch/mpx469/tGBS_enset_project/blacklists
```


### identify loci that map to contaminant sequeces in Bedadeti assembly

```
mkdir /data/scratch/mpx469/tGBS_enset_project/blacklists/blobtools

# blobtools_contigs_to_filter.txt imported to aporcrita
# contigs idenfied using blobtoolkit viewer

cut -f 6 blobtools_contigs_to_filter.txt  | grep id -v | grep -f - /data/scratch/mpx469/tGBS_enset_project/populations/populations_*_output/populations.snps.vcf

# no loci map to these regions
```


### identify loci that show high sequence homology to organlle genomes


```
mkdir /data/scratch/mpx469/tGBS_enset_project/blacklists/blastn
cd /data/scratch/mpx469/tGBS_enset_project/blacklists/blastn


# get refseq sequences

# download refseq data for chloroplast and mitochondrial sequences
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/ ftp_refseq_chloro/
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/ ftp_refseq_mito/

# gunzip
gunzip ftp_refseq_chloro/plastid.*.1.genomic.fna.gz
gunzip ftp_refseq_mito/mitochondrion.*.1.genomic.fna.gz

# spaces removed so descriptions included in output
# commas emoved just to make it look tidy
cat ftp_refseq_chloro/plastid.*.1.genomic.fna |     sed -e 's/ /_/g' -e 's/,//g' > ftp_refseq_chloro/plastid.genomic.fna
cat ftp_refseq_mito/mitochondrion.*.1.genomic.fna | sed -e 's/ /_/g' -e 's/,//g' > ftp_refseq_mito/mitochondrion.genomic.fna



# get novoplasty assembly for Bedadeti SRA data

mkdir /data/scratch/mpx469/tGBS_enset_project/blacklists/blastn/novoplasty

# cp aacross novoplasty assembled chloroplast for Bedadeti
cp /data/scratch/mpx469/sra_enset_project/novoplasty/Option_1_Bedadeti1_chloro.fasta novoplasty/

# change name of fasta sequence "Contig1" to something more meaningful
sed -i -e 's/Contig1/Bedadeti_chloro_novoplasty/g' novoplasty/Option_1_Bedadeti1_chloro.fasta



# get ensete partial chloroplast asssembly

mkdir /data/scratch/mpx469/tGBS_enset_project/blacklists/blastn/ensete

# manually download partial chloroplast assembly for ensete from:
# https://www.ncbi.nlm.nih.gov/nuccore/MH603417.1
# call "MH603417.1_Ensete_ventricosum_chloro_partial.fasta" and add to ensete folder

# spaces removed so descriptions included in output
# commas emoved just to make it look tidy
sed -i -e 's/ /_/g' -e 's/,//g' ensete/MH603417.1_Ensete_ventricosum_chloro_partial.fasta 


## get musa mitochondrial contigs

mkdir /data/scratch/mpx469/tGBS_enset_project/blacklists/blastn/musa

# download musa assembly
wget https://banana-genome-hub.southgreen.fr/sites/banana-genome-hub.southgreen.fr/files/data/fasta/version2/musa_acuminata_v2_pseudochromosome.fna -P musa/

# get mitocondrial contigs
grep -e "^>mito" musa/musa_acuminata_v2_pseudochromosome.fna | sed 's/>//g' > musa/mito.contigs

# extract mitchondrial contigs
module load seqtk
seqtk subseq musa/musa_acuminata_v2_pseudochromosome.fna musa/mito.contigs > musa/mito.contigs.fastaca

# change names to something more meaningful
sed -i -e 's/mito/musa_acuminata_v2_mito/g' musa/mito.contigs.fasta



# create blast db

cat ftp_refseq_chloro/plastid.genomic.fna \
    ftp_refseq_mito/mitochondrion.genomic.fna \
	novoplasty/Option_1_Bedadeti1_chloro.fasta \
	ensete/MH603417.1_Ensete_ventricosum_chloro_partial.fasta \
	musa/mito.contigs.fasta > organelle_seq.fasta


qsub script_makeblastdb.sh

for l in `seq 70 10 120`; do 
    for d in all_snps single_snp; do
	   echo ${l} ${d}
	done
done >> input_args

qsub script_blastn.sh

# write blast blacklists
for l in `seq 70 10 120`; do 
    for d in all_snps single_snp; do
	   cut -f 1 blast_out_${l}_${d}_top_hits | sed 's/CLocus_//g' > blacklist_blast_${l}_${d}
	done
done

```



### identify loci that show consistently high site depth

```
mkdir /data/scratch/mpx469/tGBS_enset_project/blacklists/site_depth
cd /data/scratch/mpx469/tGBS_enset_project/blacklists/site_depth



# sort and index vcf 

mkdir /data/scratch/mpx469/tGBS_enset_project/blacklists/site_depth/vcf_sorted_indexed

for l in `seq 70 10 120`; do 
    for d in all_snps single_snp; do
	   echo ${l} ${d}
	done
done >> input_args

qsub script_vcf_sort_index_array.sh



# get site depth per loci across samples

for l in `seq 70 10 120`; do 
    for d in all_snps single_snp; do
	    cat /data/scratch/mpx469/tGBS_enset_project/gstacks/popmap.txt | cut -f 1 | sed 's/.unique.sorted//g' | while read i; do
	       echo ${l} ${d} ${i}
	    done
	done
done >> input_args_site_depth

qsub script_site_depth_array.sh



# get loci info (loci chr pos) for each dataset
# required to join chr pos to loci numbers
for l in `seq 70 10 120`; do 
    for d in all_snps single_snp; do
	   grep -e "^#" -v /data/scratch/mpx469/tGBS_enset_project/populations/populations_${l}_${d}_output/populations.sumstats.tsv | cut -f 1-3 | uniq > loci_info_${l}_${d}
	done
done

grep -e "^#" -v /data/scratch/mpx469/tGBS_enset_project/populations/populations_70_all_snps_output/populations.sumstats.tsv | cut -f 1-3 > loci_info_70_all_snps


# run script to identify site with consistently high depth to be blacklisted
qsub script_identify_sites_with_high_depth_array.sh
```



### identify duplicate loci with differing 5' start sites/banana-genome-hub

```
mkdir /data/scratch/mpx469/tGBS_enset_project/blacklists/duplicates
cd /data/scratch/mpx469/tGBS_enset_project/blacklists/duplicates

for l in `seq 70 10 120`; do 
    for d in all_snps single_snp; do
	   echo ${l} ${d}
	done
done >> input_args

qsub script_duplicates_array.sh
```






## blast tGBS reads against a custom refseq bacterial database


### creating custom blast db from refseq genomes

adapted from https://github.com/dmnfarrell/epitopepredict/wiki/Creating-a-local-refseq-blast-db

```
# set up main dir for the analysis
mkdir /data/scratch/mpx469/tGBS_enset_project/blastn_refseq_bacteria
cd /data/scratch/mpx469/tGBS_enset_project/blastn_refseq_bacteria

# mkdir for ftp data
mkdir /data/scratch/mpx469/tGBS_enset_project/blastn_refseq_bacteria/ftp_refseq_bacteria

# get bacteria refseq assembly summary (accessed 10/06/2020)
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -P ftp_refseq_bacteria

# get summary information
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt -P ftp_refseq_bacteria

# get ftp files for bacterial genomes that are:
#    refererence or representative
#    Complete
awk -F "\t" '$5 == "reference genome" || $5 == "representative genome" && $12 == "Complete Genome" { print $20 } ' ftp_refseq_bacteria/assembly_summary.txt > ftp_refseq_bacteria/assembly_summary_select

wc -l ftp_refseq_bacteria/assembly_summary_select
# 3000 ftp_refseq_bacteria/assembly_summary_select

# append "*_genomic.fna.gz" to the ftp directory names
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftp_refseq_bacteria/assembly_summary_select > ftp_refseq_bacteria/assembly_summary_select_paths

# add genomes for Xcm and Xvv
echo ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/277/875/GCF_000277875.1_Xcm2005v1/GCF_000277875.1_Xcm2005v1_genomic.fna.gz >> ftp_refseq_bacteria/assembly_summary_select_paths
echo ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/015/715/GCF_003015715.1_ASM301571v1/GCF_003015715.1_ASM301571v1_genomic.fna.gz >> ftp_refseq_bacteria/assembly_summary_select_paths

qsub script_ftp_refseq_bacteria.sh
```


### blast all samples against custom refseq reference

```
# create blast databse
qsub script_makeblastdb.sh

# require taxdb in same dir to get ssciname in blast output
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar xvzf taxdb.tar.gz

# get and format ncbi lineage info

# wget new taxdump
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip

# unzip
mkdir new_taxdump && unzip new_taxdump.zip -d new_taxdump

# rm zip
rm new_taxdump.zip

# wget readme file
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/taxdump_readme.txt -P new_taxdump

# we will use "new_taxdump/rankedlineage.dmp"
# change file delimiters to make it readable in R
sed -e 's/\t|\t/|/g' -e 's/\t|//g' new_taxdump/rankedlineage.dmp > new_taxdump/rankedlineage_manual_edit.dmp

# run blast search and get read taxonomy for all samples
qsub script_blastn_taxonomy_array.sh

# create summary tables with counts per species and genera
Rscript create_summary_tables_species_genera.R
```



### count the Xanthomonas campestris pv. musacearum genome contigs/scaffolds hit in each sample 

```
cat /data/scratch/mpx469/tGBS_enset_project/sample_list.txt | while read i; do 
   echo ${i} $(grep "Xanthomonas campestris pv. musacearum NCPPB 2005" blastn_output/${i}_blastn_out | cut -f 2 | sort | uniq | wc -l)
done > summary_table_output/summary_xcm_contigs.txt 
```



### calculate proportion of the Xanthomonas campestris pv. musacearum genome identified in the blast search

```
# get ref names for all contigs in xcm assembly
grep "Xanthomonas campestris pv. musacearum" ftp_refseq_bacteria/refseq_bacteria.fasta | sed 's/>//g' | cut -f 1 -d " " > contigs_xcm

# extract Xanthomonas campestris pv. musacearum sequences
module load seqtk
sseqtk subseq ftp_refseq_bacteria/refseq_bacteria.fasta contigs_xcm | cut -f 1 -d " " > contigs_xcm.fasta

# get sequence lengths
# copied from https://www.danielecook.com/generate-fasta-sequence-lengths/
cat contigs_xcm.fasta | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > contigs_xcm_lengths

# get base pair coverage acroos all xcm contigs for each sample
qsub script_bp_coverage_per_chr_array.sh

# create a summary tables
create_summary_tables_xcm_bp_cov.R
```
