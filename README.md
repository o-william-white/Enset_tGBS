# Enset tGBS
Enset tGBS methodology

This readme file details the methodolgy used in the analysis of ensete tGBS data. All work was completed on the QMUL apocrita cluster in the follwowing directory:
```
/data/scratch/mpx469/tGBS_enset_project
```
All scripts located in the following dir
```
/data/scratch/mpx469/tGBS_enset_project/scripts
```

## Table of contents

[Import data from local hard drive into apocrita](#import-data-from-local-hard-drive-into-apocrita)

[Pre-processing of tGBS data](#pre-processing-of-tgbs-data)
   - [Create sample list to iterate through](#create-sample-list-to-iterate-through)
   - [Use trimmomatic to filter raw reads](#use-trimmomatic-to-filter-raw-reads)
   - [Use cutadapt to filter reads without restriction enzyme cutsites](#use-cutadapt-to-filter-reads-without-restriction-enzyme-cutsites)
   - [Use process radtags to truncate reads to uniform length](#use-process-radtags-to-truncate-reads-to-uniform-length)
 
 [STACKs reference mapped pipeline](#stacks-reference-mapped-pipeline)
   - [Get Bedadeti reference](#get-Bedadeti-reference)
   - [Map reads against the Bedadeti reference genome assembly](#map-reads-against-the-Bedadeti-reference-genome-assembly)
   - [Read summary statistics](#read-summary-statistics)
   - [Assemble loci with gstacks](#assemble-loci-with-gstacks)
   - [Call snps with populations](#call-snps-with-populations)

[Blacklist contaminants and paralogues](#blacklist-contaminants-and-paralogues)
   - [Identify loci that map to contaminant sequeces in Bedadeti assembly](#identify-loci-that-map-to-contaminant-sequeces-in-bedadeti-assembly)
   - [Identify loci that show high sequence homology to organelle genomes](#identify-loci-that-show-high-sequence-homology-to-organelle-genomes)
   - [Identify loci that show consistently high site depth](#identify-loci-that-show-consistently-high-site-depth)
   - [Identify duplicate loci with differing start sites](#identify-duplicate-loci-with-differing-start-sites)
   - [Create overall blacklists](#create-overall-blacklists)
   - [Repeat populations with blacklists](#repeat-populations-with-blacklists)

[Convert file formats](#convert-file-formats)

[Main phylogenetic analyses](#main-phylogenetic-analyses)
   - [Principal components analysis](#principal-components-analysis)
   - [Modeltest](#modeltest)
   - [Raxml](#raxml)
   - [Iqtree](#iqtree)
   - [Radpainter](#radpainter)
   - [Dsuite](#dsuite)
   - [EasySFS](#easysfs)

[Repeat SNP calling for population genetic dataset](#repeat-snp-calling-for-population-genetic-dataset)



<br/>
<div align="right">
    <b><a href="#enset-tgbs">↥ back to top</a></b>
</div>
<br/>

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
    <b><a href="#enset-tgbs">↥ back to top</a></b>
</div>
<br/>

## Pre processing of tGBS data


### Create sample list to iterate through

```
# set dir
cd /data/scratch/mpx469/tGBS_enset_project
```

Import GBS_metadata in .csv format and create a sample list to iterate through
```
ls -1 Data2Bio_final/raw/*fq.gz | sed -e 's/.digested.fq.gz//g' -e 's:Data2Bio_final/raw/::g' > sample_list.txt
wc -l sample_list.txt
# 283 sample_list.txt
```



### Use trimmomatic to filter raw reads

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/trimmomatic
cd /data/scratch/mpx469/tGBS_enset_project/trimmomatic

# cp script
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_trimmomatic_array.sh .

# run trimmomatic
qsub script_trimmomatic_array.sh
```



### Use cutadapt to filter reads without restriction enzyme cutsites

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/cutadapt
cd /data/scratch/mpx469/tGBS_enset_project/cutadapt

# cp script
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_cutadapt_array.sh .

# run cutadapt
qsub script_cutadapt_array.sh
```



### Use process radtags to truncate reads to uniform length

Will truncate reads at varying read lengths (60 bp to 120 bp in steps of 10)

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/process_radtags
cd /data/scratch/mpx469/tGBS_enset_project/process_radtags

# create input args for array script
for l in `seq 70 10 120`; do 
   cat /data/scratch/mpx469/tGBS_enset_project/sample_list.txt | while read i; do
      echo $l $i >> input_args
   done	  
done

wc -l input_args
# 1698 input_args

# cp script
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_process_radtags_array.sh .

# run process radtags
qsub script_process_radtags_array.sh
```




<br/>
<div align="right">
    <b><a href="#enset-tgbs">↥ back to top</a></b>
</div>
<br/>

## STACKs reference mapped pipeline

### Get Bedadeti reference
```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/bedadeti_assembly
cd /data/scratch/mpx469/tGBS_enset_project/bedadeti_assembly

# cp script
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_wget_Bedadeti.sh .

# run script
qsub script_wget_Bedadeti.sh
```


### Map reads against the Bedadeti reference genome assembly

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti
cd /data/scratch/mpx469/tGBS_enset_project/map_reads_bedadeti

# cp scripts
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_bwa_index.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_map_reads_bedadeti_array.sh .

# cp input_args to dir
cp /data/scratch/mpx469/tGBS_enset_project/process_radtags/input_args .

# index genome 
qsub script_bwa_index.sh

# map reads
qsub script_map_reads_bedadeti_array.sh
```



### Read summary statistics
```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/summary_stats/
cd /data/scratch/mpx469/tGBS_enset_project/summary_stats/

# cp script
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_read_summary_array.sh .

# set up output prior to running script
for LENGTH in `seq 70 10 120`; do
   echo SAMPLE,RAW,TRIMMOMATIC,CUTADAPT,CUTADAPT_BP,PROCESS_RADTAGS,PROCESS_RADTAGS_BP,BWA_UNMAPPED,BWA_MAPPED,BWA_XA_SA,SAM_UNMAPPED,SAM_MAPPED > summary_$LENGTH.csv
done

# run script
qsub script_read_summary_array.sh
```



### Assemble loci with gstacks 

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/gstacks
cd /data/scratch/mpx469/tGBS_enset_project/gstacks

# cp scripts
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_gstacks_array.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/write_popmap.R .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/plot_gstacks_summary.R .

# write popmap
# note that this requires the .csv file in /data/scratch/mpx469/tGBS_enset_project/tGBS_metadata_phylogenetic_analysis.csv
module add R
Rscript write_popmap.R

# run gstacks
qsub script_gstacks_array.sh

# get summary stats for loci assembled at each read length
echo -e 'loci reads coverage' > summary_gstacks

for l in `seq 70 10 120`; do 
   echo ${l} $(grep Built gstacks_${l}_output/gstacks.log | cut -f 2,5 -d " ") $(grep coverage gstacks_${l}_output/gstacks.log | cut -f 6 -d " " | sed -e 's/mean=//g' -e 's/x,//g' )  >> summary_gstacks
done

# plot summary stats
Rscript plot_gstacks_summary.R

```



### Call snps with populations

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/populations
cd /data/scratch/mpx469/tGBS_enset_project/populations

# cp scripts
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_populations_all_snps_array.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_populations_single_snp_array.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/plot_populations_summary.R .

# run populations
qsub script_populations_all_snps_array.sh
qsub script_populations_single_snp_array.sh

# get summary stats for each assembly
echo -e 'loci sites filtered variant' > summary_all_snps
echo -e 'loci sites filtered variant' > summary_single_snp

for l in `seq 70 10 120`; do 
    echo ${l} $(grep Kept populations_${l}_all_snps_output/populations.log | cut -f 2,6,8,14 -d " ") >> summary_all_snps
    echo ${l} $(grep Kept populations_${l}_single_snp_output/populations.log | cut -f 2,6,8,14 -d " ") >> summary_single_snp
done

# summary plot
Rscript plot_populations_summary.R

```


<br/>
<div align="right">
    <b><a href="#enset-tgbs">↥ back to top</a></b>
</div>
<br/>

## Blacklist contaminants and paralogues

```
mkdir /data/scratch/mpx469/tGBS_enset_project/blacklists
```

### Identify loci that map to contaminant sequeces in Bedadeti assembly

```
mkdir /data/scratch/mpx469/tGBS_enset_project/blacklists/blobtools
cd /data/scratch/mpx469/tGBS_enset_project/blacklists/blobtools

# blobtools_contigs_to_filter.txt imported to apocrita
# contigs identfied using blobtoolkit viewer

# get contigs 
cut -f 6 blobtools_contigs_to_filter.txt  | grep id -v > contigs.txt

# grep for contigs in vcf files
for l in `seq 70 10 120`; do 
    for d in all_snps single_snp; do
	echo Checking vcf for read length $i and dataset ${d}
   grep -f contigs.txt /data/scratch/mpx469/tGBS_enset_project/populations/populations_${l}_${d}_output/populations.snps.vcf   
    done
done

# no loci map to these regions
```


### Identify loci that show high sequence homology to organelle genomes


```
mkdir /data/scratch/mpx469/tGBS_enset_project/blacklists/blastn
cd /data/scratch/mpx469/tGBS_enset_project/blacklists/blastn

# set up blast database
mkdir /data/scratch/mpx469/tGBS_enset_project/blacklists/blastn/references


# get refseq sequences

# download refseq data for chloroplast and mitochondrial sequences
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/       ftp_refseq_chloro/
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/ ftp_refseq_mito/

# gunzip
gunzip ftp_refseq_chloro/plastid.*.1.genomic.fna.gz
gunzip ftp_refseq_mito/mitochondrion.*.1.genomic.fna.gz

# cat ftp download to a single file
#	spaces removed so descriptions included in output
#	commas emoved just to make it look tidy
#	there was one assembly with a "#" in name which was also removed
cat ftp_refseq_chloro/plastid.*.1.genomic.fna     | sed -e 's/ /_/g' -e 's/,//g' -e 's/#/_/' > references/refseq_plastid_genomic.fasta
cat ftp_refseq_mito/mitochondrion.*.1.genomic.fna | sed -e 's/ /_/g' -e 's/,//g' -e 's/#/_/' > references/refseq_mitochondrion_genomic.fasta

# rm ftp download
rm -r ftp_refseq_chloro
rm -r ftp_refseq_mito


# get novoplasty assembly for Bedadeti SRA data

# cp across novoplasty assembled chloroplast for Bedadeti
cp /data/scratch/mpx469/assemble_bedadeti_plastome/novoplasty/Option_1_Bedadeti1.fasta references/

# change name of fasta sequence "Contig1" to something more meaningful
sed -i -e 's/Contig1/Bedadeti_chloro_novoplasty/g' references/Option_1_Bedadeti1.fasta


# get ensete partial chloroplast asssembly

# manually download partial chloroplast assembly for ensete from:
# https://www.ncbi.nlm.nih.gov/nuccore/MH603417.1
# call "MH603417.1_Ensete_ventricosum_chloro_partial.fasta" and add to references folder

# spaces removed so descriptions included in output
# commas emoved just to make it look tidy
sed -i -e 's/ /_/g' -e 's/,//g' references/MH603417.1_Ensete_ventricosum_chloro_partial.fasta 


# get musa mitochondrial contigs

# download musa assembly
# note this webpage no longer available # replaced by https://banana-genome-hub.southgreen.fr/node/50/413 and I was unable to find the relevant file
wget https://banana-genome-hub.southgreen.fr/sites/banana-genome-hub.southgreen.fr/files/data/fasta/version2/musa_acuminata_v2_pseudochromosome.fna -P references/

# get mitocondrial contigs
grep -e "^>mito" references/musa_acuminata_v2_pseudochromosome.fna | sed 's/>//g' > references/mito_contigs

# extract mitchondrial contigs
module load seqtk
seqtk subseq references/musa_acuminata_v2_pseudochromosome.fna references/mito_contigs > references/musa_acuminata_mito_contigs.fasta 

# change names to something more meaningful
sed -i -e 's/mito/musa_acuminata_v2_mito/g' references/musa_acuminata_mito_contigs.fasta


# create blast db
cat references/*.fasta > organelle.fasta

# cp scripts
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_makeblastdb_organelle.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_blastn_organelle.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/top_hit.R .

# makeblastdb
qsub script_makeblastdb_organelle.sh

# create inupt list
for l in `seq 70 10 120`; do 
    for d in all_snps single_snp; do
	   echo ${l} ${d}
	done
done >> input_args

# run blastn and write blacklists
qsub script_blastn.sh
```



### Identify loci that show consistently high site depth

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/blacklists/loci_depth
cd /data/scratch/mpx469/tGBS_enset_project/blacklists/loci_depth

# cp script
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_identify_high_depth_loci.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/identify_high_depth_loci.R .

# create inupt list
for l in `seq 70 10 120`; do 
    for d in all_snps single_snp; do
	   echo ${l} ${d}
	done
done >> input_args

# run script
qsub script_identify_high_depth_loci.sh
```



### Identify duplicate loci with differing start sites

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/blacklists/duplicate_loci
cd /data/scratch/mpx469/tGBS_enset_project/blacklists/duplicate_loci

# create input list
for l in `seq 70 10 120`; do 
    for d in all_snps single_snp; do
	   echo ${l} ${d}
	done
done >> input_args

# cp script
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_identify_duplicate_loci.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/identify_duplicate_loci.R .

# run script
qsub script_identify_duplicate_loci.sh
```


### Create overall blacklists

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/blacklists/blacklists_overall
cd /data/scratch/mpx469/tGBS_enset_project/blacklists/blacklists_overall

# cp script
cp /data/scratch/mpx469/tGBS_enset_project/scripts/write_overall_blacklist.R .

# run Rscript
Rscript write_overall_blacklist.R
```


### Repeat populations with blacklists

```
# set dir
cd /data/scratch/mpx469/tGBS_enset_project/populations

# cp scripts
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_populations_all_snps_blacklist_array.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_populations_single_snp_blacklist_array.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/plot_populations_blacklist_summary.R .

# run populations
qsub script_populations_all_snps_blacklist_array.sh
qsub script_populations_single_snp_blacklist_array.sh

# get summary stats for each assembly
echo -e 'loci sites filtered variant' > summary_all_snps_blacklist
echo -e 'loci sites filtered variant' > summary_single_snp_blacklist

for l in `seq 70 10 120`; do 
    echo ${l} $(grep Kept populations_${l}_all_snps_blacklist_output/populations.log | cut -f 2,6,8,14 -d " ") >> summary_all_snps_blacklist
    echo ${l} $(grep Kept populations_${l}_single_snp_blacklist_output/populations.log | cut -f 2,6,8,14 -d " ") >> summary_single_snp_blacklist
done

# summary plot
Rscript plot_populations_blacklist_summary.R
```

<br/>
<div align="right">
    <b><a href="#enset-tgbs">↥ back to top</a></b>
</div>
<br/>



## Convert file formats

The stacks phylip output is in interleaved format. For raxml-ng it need to be in sequential format. We also create other potentailly useful formats including nexus and fasta.

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/convert_file_formats
cd /data/scratch/mpx469/tGBS_enset_project/convert_file_formats

# cp scripts
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_convert_file_formats_array.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/convert_alignment_format.py .

# create input list
for l in `seq 70 10 120`; do 
    for d in all_snps single_snp; do
	   echo ${l} ${d}
	done
done >> input_args

qsub script_convert_file_formats_array.sh
```

<br/>
<div align="right">
    <b><a href="#enset-tgbs">↥ back to top</a></b>
</div>
<br/>


## Main phylogenetic analyses

### Principal components analysis

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/pca
cd /data/scratch/mpx469/tGBS_enset_project/pca

# cp scripts
cp /data/scratch/mpx469/tGBS_enset_project/scripts/plot_plink_pca.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/plot_plink_pca.R .

# need a tab-delimited text file with family IDs in the first column and within-family IDs for outgroup samples
grep Outgroup /data/scratch/mpx469/tGBS_enset_project/tGBS_metadata_phylogenetic_analysis.csv | cut -f2 -d "," | grep -f - /data/scratch/mpx469/tGBS_enset_project/populations/populations_80_single_snp_blacklist_output/populations.plink.ped | cut -f 1,2 > outgroups.txt

# run script
bash plot_plink_pca.sh
```




### Modeltest

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/modeltest_ng
cd /data/scratch/mpx469/tGBS_enset_project/modeltest_ng

# cp script
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_modeltest_ng_array.sh .

# run modeltest
qsub script_modeltest_ng_array.sh
```




### Raxml

Maximum likelihood tree estimation with raxml-ng. We analyse the dataset with a read length of 80 and a single snp per locus.

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/raxml_ng
mkdir /data/scratch/mpx469/tGBS_enset_project/raxml_ng/raxml_ng_80_single_snp
cd /data/scratch/mpx469/tGBS_enset_project/raxml_ng/raxml_ng_80_single_snp

# cp scripts
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_raxml_ng_* .

# run script to parse alignment
qsub script_raxml_ng_01_parse.sh

# run scripts for tree searches
qsub script_raxml_ng_02_tree_search_rand_array.sh
qsub script_raxml_ng_03_tree_search_pars_array.sh
qsub script_raxml_ng_04_bootstrap_array.sh

# check if jobs completed
cat raxml_ng_output/all_tree_search_rand*.log | grep -e "Elapsed time" -c
cat raxml_ng_output/all_tree_search_pars*.log | grep -e "Elapsed time" -c
cat raxml_ng_output/all_bootstrap*.log | grep -e "Elapsed time" -c

# resubmit jobs if they did not complete
# raxml-ng will retart from checkpoint
# note these jobs are submitted individually rather than as an array

# rand tree search
for i in {1..5000}; do
   RAND=$(grep "Elapsed time" raxml_ng_output/all_tree_search_rand_${i}.raxml.log -c)
   if [ $RAND -eq 0 ]; then
      qsub -t ${i} script_raxml_ng_02_tree_search_rand_array.sh
   fi
done

# pars tree search
for i in {1..5000}; do
   PARS=$(grep "Elapsed time" raxml_ng_output/all_tree_search_pars_${i}.raxml.log -c)
   if [ $PARS -eq 0 ]; then
      qsub -t ${i} script_raxml_ng_03_tree_search_pars_array.sh
   fi
done

# bootstrap
for i in {1..5000}; do
   BOOT=$(grep "Elapsed time" raxml_ng_output/all_bootstrap_${i}.raxml.log -c)
   if [ $BOOT -eq 0 ]; then
      qsub -t ${i} script_raxml_ng_04_bootstrap_array.sh
   fi
done

# check again that all jobs have completed as above

# see best scoring trees
grep "Final LogLikelihood:" raxml_ng_output/all_tree_search_*.raxml.log | sort -k 3 | head

# create symbolic link to best tree
ln -s `grep "Final LogLikelihood:" raxml_ng_output/all_tree_search_*.raxml.log | sort -k 3 | head -n 1 | cut -f 1 -d ":" | sed 's/log/bestTree/g'` best.tre

# cat all bootstrap trees
cat raxml_ng_output/all_bootstrap*bootstraps > bootstrap_trees

# check bootstrap convergence
qsub script_raxml_ng_05_bsconverge.sh

# add bootstrap values to best tree

module add gcc/7.1.0
export PATH=/data/home/mpx469/software/raxml-ng-pthreads/raxml-ng/bin/:$PATH

raxml-ng \
   --support \
   --tree best.tre \
   --bs-trees bootstrap_trees \
   --bs-metric fbp \
   --prefix best.tre.fbp

raxml-ng \
   --support \
   --tree best.tre \
   --bs-trees bootstrap_trees \
   --bs-metric tbe \
   --prefix best.tre.tbe

# root best supported trees
cp /data/scratch/mpx469/tGBS_enset_project/scripts/root_and_ladderise_tree.py .
source /data/home/mpx469/software/python-virtualenv/ete3/bin/activate
python root_and_ladderise_tree.py best.tre.fbp.raxml.support best.tre.fbp.raxml.support.rooted.newick pop82,pop160,pop162
python root_and_ladderise_tree.py best.tre.tbe.raxml.support best.tre.tbe.raxml.support.rooted.newick pop82,pop160,pop162
```




### Iqtree

Maximum likelihood tree estimation with iq-tree. We use the same data as above for raxml-ng. 

```
# set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/iq_tree
mkdir /data/scratch/mpx469/tGBS_enset_project/iq_tree/iq_tree_80_single_snp
cd /data/scratch/mpx469/tGBS_enset_project/iq_tree/iq_tree_80_single_snp

# cp scripts
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_iq_tree_array.sh .

# write input args
for i in `seq 0.1 0.1 0.5`; do
   for r in {1..100}; 
      do echo -pers ${i} -nstop 1000 -pre iq_tree_output/out_pers_${i}_r${r}; 
   done 
done > input_args

# run iqtree
qsub script_iq_tree_array.sh

# note I needed to repeat some scripts, but it is no problem with checkpoints

# get tree likelihood scores
grep -e "BEST SCORE FOUND" iq_tree_output/*log | sed -e 's/:BEST SCORE FOUND : /\t/g' -e 's/.log/.treefile/g' > tree_likelihood.txt

# identify tree with the best score
BEST_TREE=$(sort -k 2 -n -r tree_likelihood.txt | head -n 1 | cut -f 1)

# create symbolic link
ln -s $BEST_TREE best.tre.iqtree
```
Tree copied to local dir and rooted on outgroups pop82,pop160,pop162 and ladderised




### Radpainter

Note for this analysis we use a clone corrected dataset, with representatives of each MLG in the file "mlg_farthest_bitwise_monophyletic_single_rep.txt"

```
#set up dir
mkdir /data/scratch/mpx469/tGBS_enset_project/radpainter
cd    /data/scratch/mpx469/tGBS_enset_project/radpainter

# cp scripts
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_radpainter.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/clone_correct_radpainter_input.R . 
cp /data/scratch/mpx469/tGBS_enset_project/scripts/mlg_farthest_bitwise_monophyletic_single_rep.txt . 

# run radpainter
qsub script_radpainter.sh
```




### Dsuite
```
#set up dir
mkdir /data/scratch/mpx469/tGBS_enset_project/dsuite
cd    /data/scratch/mpx469/tGBS_enset_project/dsuite

# cp scripts
cp /data/scratch/mpx469/tGBS_enset_project/scripts/write_dsuite_sets.R .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/write_tree_renamed_tip_labels.R .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_dsuite_by_pop.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_dsuite_by_sample.sh .

# mlg_farthest_bitwise_monophyletic_single_rep.txt copied from local dir

# run scripts
Rscript write_dsuite_sets.R
Rscript write_tree_renamed_tip_labels.R
qsub script_dsuite_by_pop.sh
qsub script_dsuite_by_sample.sh
```




### EasySFS
Set up data for fastsimcoal
```
mkdir /data/scratch/mpx469/tGBS_enset_project/easy_sfs
cd    /data/scratch/mpx469/tGBS_enset_project/easy_sfs

module load anaconda3/
conda activate easySFS
export PATH=/data/home/mpx469/software/conda/easySFS/:$PATH

# popmap_sfs.txt copied from local treemix dir
# the popmap contain samples data for 1 wild and 5 domesticated populations
# select 1 wild and 2 domesticated
grep -e Wild -e Domesticated_2 -e Domesticated_5 popmap_sfs.txt > popmap_sfs_select.txt

# cp vcf to dir
cp /data/scratch/mpx469/tGBS_enset_project/populations/populations_80_single_snp_blacklist_output/populations.snps.vcf .

# we cannot use a folded sfs as we do not have accurately polarised data 
# i.e. we do not know the ancestral allele

# preview
easySFS.py -i populations.snps.vcf -p popmap_sfs_select.txt --order Wild,Domesticated_2,Domesticated_5 --preview

# Continue, excluding samples not in both pops file and VCF? (yes/no)
# yes

#Wild
#(2, 1739)       (3, 2609)       (4, 3192)       (5, 3625)       (6, 3975)       (7, 4252)       (8, 4493)       (9, 4687)       (10, 4864)      (11, 4980)      (12, 5115)      (13, 5164)      (14, 5268)      (15, 5232)      (16, 5313)(17, 5131)      (18, 5196)      (19, 4802)      (20, 4853)      (21, 4055)      (22, 4091)      (23, 2896)      (24, 2918)      (25, 1314)      (26, 1323)
#
#
#Domesticated_2
#(2, 1087)       (3, 1630)       (4, 1984)       (5, 2243)       (6, 2444)       (7, 2608)       (8, 2746)       (9, 2863)       (10, 2965)      (11, 3054)      (12, 3134)      (13, 3197)      (14, 3262)      (15, 3276)      (16, 3331)(17, 3272)      (18, 3318)      (19, 3038)      (20, 3075)      (21, 2476)      (22, 2503)      (23, 1617)      (24, 1635)      (25, 694)       (26, 702)
#
#
#Domesticated_5
#(2, 794)        (3, 1192)       (4, 1444)       (5, 1623)       (6, 1760)       (7, 1870)       (8, 1961)       (9, 2038)       (10, 2105)      (11, 2160)      (12, 2213)      (13, 2260)      (14, 2303)      (15, 2332)      (16, 2369)(17, 2370)      (18, 2402)      (19, 2322)      (20, 2348)      (21, 2104)      (22, 2126)      (23, 1540)      (24, 1555)      (25, 765)       (26, 773)
#

Each column is the number of samples in the projection and the number of segregating sites at that projection value. 
The dadi manual recommends maximizing the number of segregating sites, but at the same time if you have lots of missing
data then you might have to balance # of segregating sites against # of samples to avoid downsampling too far.

# max segregating sites
# Wild           (16, 5313)
# Domesticated_2 (16, 3331)
# Domesticated_5 (18, 2399)

# easysfs
easySFS.py -i populations.snps.vcf -p popmap_sfs_select.txt --order Wild,Domesticated_2,Domesticated_5 --proj 16,16,18

# Continue, excluding samples not in both pops file and VCF? (yes/no)
# yes

# note pop order is very important 
# 0 = Wild
# 1 = Domesticated_2
# 2 = Domesticated_5
```

<br/>
<div align="right">
    <b><a href="#enset-tgbs">↥ back to top</a></b>
</div>
<br/>


### Repeat SNP calling for population genetic dataset

See README_population_genetic_dataset.md for details
