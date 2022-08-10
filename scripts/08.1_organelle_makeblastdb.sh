#!/bin/sh
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -N job_organelle_make_blastdb

ml entrez
ml seqtk
ml blast+

mkdir -p blacklist
mkdir -p blacklist/blastn
mkdir -p blacklist/blastn/references

# download refseq data for chloroplast and mitochondrial sequences
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/       blacklist/blastn/references/ftp_refseq_chloro/
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/ blacklist/blastn/references/ftp_refseq_mito/

# gunzip
gunzip blacklist/blastn/references/ftp_refseq_chloro/plastid.*.1.genomic.fna.gz
gunzip blacklist/blastn/references/ftp_refseq_mito/mitochondrion.*.1.genomic.fna.gz

# cat ftp download to a single file
#   spaces removed so descriptions included in output
#   commas emoved just to make it look tidy
#   there was one assembly with a "#" in name which was also removed
cat blacklist/blastn/references/ftp_refseq_chloro/plastid.*.1.genomic.fna     | sed -e 's/ /_/g' -e 's/,//g' -e 's/#/_/' > blacklist/blastn/references/refseq_plastid_genomic.fasta
cat blacklist/blastn/references/ftp_refseq_mito/mitochondrion.*.1.genomic.fna | sed -e 's/ /_/g' -e 's/,//g' -e 's/#/_/' > blacklist/blastn/references/refseq_mitochondrion_genomic.fasta

# rm ftp download
rm -r blacklist/blastn/references/ftp_refseq_chloro
rm -r blacklist/blastn/references/ftp_refseq_mito

# get novoplasty assembly for Bedadeti SRA data
cp additional_data/Option_1_Bedadeti1.fasta blacklist/blastn/references

# get ensete partial chloroplast asssembly
efetch -db nucleotide -id MH603417.1 -format fasta > blacklist/blastn/references/MH603417.1_Ensete_ventricosum_chloro_partial.fasta

# get musa chloroplast assemblies
efetch -db nucleotide -id NC_042874.1 -format fasta > blacklist/blastn/references/NC_042874.1_Musa_ornata_chloro.fasta
efetch -db nucleotide -id NC_039815.1 -format fasta > blacklist/blastn/references/NC_039815.1_Musa_balbisiana_var_balbisiana_chloro.fasta
efetch -db nucleotide -id NC_022926.1 -format fasta > blacklist/blastn/references/NC_022926.1_Musa_textilis_chloro.fasta

# spaces and commas removed from description
sed -i -e 's/ /_/g' -e 's/,//g' blacklist/blastn/references/*.fasta

# download musa assembly
wget https://banana-genome-hub.southgreen.fr/filebrowser/download/444/musa_acuminata_v2_pseudochromosome.fna -P blacklist/blastn/references/

# get mitocondrial contigs
grep -e "^>mito" blacklist/blastn/references/musa_acuminata_v2_pseudochromosome.fna | sed 's/>//g' > blacklist/blastn/references/mito_contigs

# extract mitchondrial contigs
seqtk subseq blacklist/blastn/references/musa_acuminata_v2_pseudochromosome.fna blacklist/blastn/references/mito_contigs > blacklist/blastn/references/musa_acuminata_mito_contigs.fasta 

# change names to something more meaningful
sed -i -e 's/mito/musa_acuminata_v2_mito/g' blacklist/blastn/references/musa_acuminata_mito_contigs.fasta

# rm intermediate
rm blacklist/blastn/references/mito_contigs
rm blacklist/blastn/references/musa_acuminata_v2_pseudochromosome.fna

# list dir contents
ls -lh blacklist/blastn/references/

# cat contents
cat blacklist/blastn/references/*.fasta > blacklist/blastn/organelle.fasta

# makebdb
makeblastdb -in blacklist/blastn/organelle.fasta -out blacklist/blastn/organelle.db -dbtype nucl

echo done
