
# Repeat snp calling for a population genetic dataset

### cp new popmaps to gstacks dir
```
# set dir
cd /data/scratch/mpx469/tGBS_enset_project/gstacks
```
created popmap with clone corrected domesticated and wild samples in 
"C:\Users\owh10kg\Desktop\tGBS_enset_project\estimate_mlg"
popmap_select* copied to dir


### run populaitons 
```
# set dir
cd /data/scratch/mpx469/tGBS_enset_project/populations

# create input args for array
for l in `seq 70 10 120`; do 
	for p in all dom wil; do
	   echo -e "${l}\t${p}"
	done
done > input_args_select

# cp scripts
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_populations_select_single_snp_array.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_populations_select_all_snps_array.sh .

# run populations
qsub script_populations_select_all_snps_array.sh
qsub script_populations_select_single_snp_array.sh
```

### write whitelist for each dataset

Note that the single snps dataset is treated differently to the all snps dataset

For the single_snp dataset, it is possible for a site to be heterozygous for one category (wil) and not for another (dom)
in this situation another snp will be selelcted for the dom and the resulting whitelist has more than one snp for a given locus

Hence we filter all_snps for uniq loci_pos combinations and filter single_snp for the first snp in the whitelist

```
# create if else statement within loop
# if all_snps take all uniq sites
# if single_snp take onlt the first snp in a locus

for l in `seq 70 10 120`; do 
    for d in all_snps single_snp; do

       if [ $d == all_snps ]

       then
          
          # create overall whitelist
          cat populations_select_${l}_${d}_output/all/whitelist \
              populations_select_${l}_${d}_output/dom/whitelist \
              populations_select_${l}_${d}_output/wil/whitelist | grep -e "Locus ID" -v | sort -n | uniq > populations_select_${l}_${d}_output/whitelist_overall

       else

          cat populations_select_${l}_${d}_output/all/whitelist \
              populations_select_${l}_${d}_output/dom/whitelist \
              populations_select_${l}_${d}_output/wil/whitelist | grep -e "Locus ID" -v | sort -k 1 -n -u > populations_select_${l}_${d}_output/whitelist_overall

       fi

   done
done

# number of loci recovered by single and all snps is the same
for l in `seq 70 10 120`; do 
   for d in all_snps single_snp; do
      echo ${l} ${d} $(cut -f 1 populations_select_${l}_${d}_output/whitelist_overall | sort -n | uniq | wc -l)
   done
done

# and there are no duplicate loci in the single snp dataset
for l in `seq 70 10 120`; do 
   echo ${l} $(cut -f 1 populations_select_${l}_single_snp_output/whitelist_overall | sort -n | uniq -d | wc -l)
done

# cp scripts
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_populations_select_single_snp_whitelist_array.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_populations_select_all_snps_whitelist_array.sh .

# run populations whitelist
qsub script_populations_select_single_snp_whitelist_array.sh
qsub script_populations_select_all_snps_whitelist_array.sh

# summary populations whitelist
echo -e 'loci sites filtered variant' > summary_select_all_snps_whitelist
echo -e 'loci sites filtered variant' > summary_select_single_snp_whitelist

for l in `seq 70 10 120`; do 
    echo ${l} $(grep Kept populations_select_${l}_all_snps_whitelist_output/populations.log | cut -f 2,6,8,14 -d " ") >> summary_select_all_snps_whitelist
    echo ${l} $(grep Kept populations_select_${l}_single_snp_whitelist_output/populations.log | cut -f 2,6,8,14 -d " ") >> summary_select_single_snp_whitelist
done
```

Compare summary to whitelist length

all_snps 
variants ~ whitelist length

```
cat summary_select_all_snps_whitelist
```
loci sites filtered variant
70 26467 1927349 4 62163
80 26596 2210616 3 69627
90 25643 2394628 6 73097
100 24188 2507480 5 75259
110 23205 2644152 5 78496
120 21894 2717029 1 80101

```
for l in `seq 70 10 120`; do 
    echo ${l} $(wc -l populations_select_${l}_all_snps_output/whitelist_overall | cut -f 1 -d " ")
done
```
70 62167
80 69630
90 73103
100 75264
110 78501
120 80102


single_snp 
loci ~ whitelist length

```
cat summary_select_single_snp_whitelist
```
70 26467 1927349 4 26463
80 26596 2210616 2 26594
90 25643 2394628 6 25637
100 24188 2507480 1 24187
110 23205 2644152 4 23201
120 21894 2717029 1 21893

```
for l in `seq 70 10 120`; do 
    echo ${l} $(wc -l populations_select_${l}_single_snp_output/whitelist_overall| cut -f 1 -d " ")
done
```

some loci filtered from those specified in whitelist

stacks website 
http://catchenlab.life.illinois.edu/stacks/manual/#pfiles

Why do loci drop out of the analysis despite being on the whitelist?
Loci, or SNPs within a locus can still drop out from an analysis despite 
being on the whitelist. 
This can happen for several reasons, including:

The filters in the populations program are still applied after the whitelist
is applied. A locus must still pass the filters to be retained. Once you 
have created a whitelist, it is normal to turn off most or all other filters.
If you change your population map after creating the whitelist, you may see
SNPs drop out of the analysis because introducing a population map may change
if a locus is fixed. In a large, single population a locus may be polymorphic,
but once you subset your data into multiple populations that locus may become 
fixed in one or more subpopulations and will not be output in those populations.
If you add populations to your population map, you may find a small number of 
loci where additional populations bring a third or fourth allele into a particular
SNP position, causing that position to fail the infinite alleles assumption 
of the software and be dropped.



### Blacklist contaminants and paralogues

#### Identify loci that map to contaminant sequeces in Bedadeti assembly
```
# set dir 
cd /data/scratch/mpx469/tGBS_enset_project/blacklists/blobtools

# grep for contigs in vcf files
cut -f 6 blobtools_contigs_to_filter.txt  | grep id -v | grep -f - /data/scratch/mpx469/tGBS_enset_project/populations/populations_select_*_whitelist_output/populations.snps.vcf
```
No loci map to these regions


#### Identify loci that show high sequence homology to organelle genomes
```
# set dir
cd /data/scratch/mpx469/tGBS_enset_project/blacklists/blastn

# cp script
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_blastn_organelle_select.sh .

# run blast
qsub script_blastn_organelle_select.sh
```

#### Identify loci that show consistently high site depth

```
# set dir
cd /data/scratch/mpx469/tGBS_enset_project/blacklists/loci_depth

# cp script
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_identify_high_depth_loci_select.sh .

# run script
qsub script_identify_high_depth_loci_select.sh
```


#### Identify duplicate loci with differing start sites

```
# set dir
cd /data/scratch/mpx469/tGBS_enset_project/blacklists/duplicate_loci

# cp script


# run script
qsub script_identify_duplicate_loci_select.sh
```

#### Create overall blacklists

```
# set dir
cd /data/scratch/mpx469/tGBS_enset_project/blacklists/blacklists_overall

# cp script
cp /data/scratch/mpx469/tGBS_enset_project/scripts/write_overall_blacklist_select.R .

# run Rscript
Rscript Rscript write_overall_blacklist_select.R
```

#### Repeat populations with blacklists
```
# set dir
cd /data/scratch/mpx469/tGBS_enset_project/populations

# cp scripts
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_populations_select_all_snps_whitelist_blacklist_array.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_populations_select_single_snp_whitelist_blacklist_array.sh .

# run populations
qsub script_populations_select_all_snps_whitelist_blacklist_array.sh
qsub script_populations_select_single_snp_whitelist_blacklist_array.sh

# summary populations whitelist
echo -e 'loci sites filtered variant' > summary_select_all_snps_whitelist_blacklist
echo -e 'loci sites filtered variant' > summary_select_single_snp_whitelist_blacklist

for l in `seq 70 10 120`; do 
    echo ${l} $(grep Kept populations_select_${l}_all_snps_whitelist_blacklist_output/populations.log   | cut -f 2,6,8,14 -d " ") >> summary_select_all_snps_whitelist_blacklist
    echo ${l} $(grep Kept populations_select_${l}_single_snp_whitelist_blacklist_output/populations.log | cut -f 2,6,8,14 -d " ") >> summary_select_single_snp_whitelist_blacklist
done
```


### Annotate SNPs with snpeff

These commands were run on an interative node 

```
set dir
mkdir /data/scratch/mpx469/tGBS_enset_project/snpeff
mkdir /data/scratch/mpx469/tGBS_enset_project/snpeff/data
mkdir /data/scratch/mpx469/tGBS_enset_project/snpeff/data/genomes
mkdir /data/scratch/mpx469/tGBS_enset_project/snpeff/data/bedadeti

# Get annotation files
cd /data/scratch/mpx469/tGBS_enset_project/snpeff/data/bedadeti
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/818/735/GCA_000818735.3_Bedadeti_annotated/GCA_000818735.3_Bedadeti_annotated_genomic.gtf.gz
mv GCA_000818735.3_Bedadeti_annotated_genomic.gtf.gz genes.gtf.gz

# Get the genome
cd /data/scratch/mpx469/tGBS_enset_project/snpeff/data/genomes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/818/735/GCA_000818735.3_Bedadeti_annotated/GCA_000818735.3_Bedadeti_annotated_genomic.fna.gz
mv GCA_000818735.3_Bedadeti_annotated_genomic.fna.gz bedadeti.fa.gz

# cp scripts from software download
cd /data/scratch/mpx469/tGBS_enset_project/snpeff
cp /data/home/mpx469/software/SnpEff/snpEff/snpEff.jar .
cp /data/home/mpx469/software/SnpEff/snpEff/SnpSift.jar .
cp /data/home/mpx469/software/SnpEff/snpEff/snpEff.config .
cp /data/home/mpx469/software/SnpEff/snpEff/scripts/vcfEffOnePerLine.pl .

# add entry to config file
echo "bedadeti.genome : Enset" >> snpEff.config 

# cp scripts
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_build_database.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_snpeff_ann.sh .
cp /data/scratch/mpx469/tGBS_enset_project/scripts/script_snpeff_eff.sh .

# cp populations input
cp /data/scratch/mpx469/tGBS_enset_project/populations/populations_select_80_all_snps_whitelist_blacklist_output/populations.snps.vcf .

# build database
qsub script_build_database.sh

# run snnotation scripts
qsub script_snpeff_ann.sh
qsub script_snpeff_eff.sh
```


