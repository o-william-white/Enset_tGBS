
# download blobtoolkit data bedadeti
# go to blobtoolkit viewer webpage https://blobtoolkit.genomehubs.org/view/
# enter "Ensete ventricosum" to search bar
# select Dataset JTFG03, accession GCA_000818735.3 with 45,742 scaffolds
# select "table" tab
# select csv download
# saved to additional_data/blobtoolkit_GCA_000818735.3.csv

mkdir -p blacklist
mkdir -p blacklist/select_contaminant

tail -n +2 additional_data/blobtoolkit_GCA_000818735.3.csv | grep -e "Streptophyta" -e "no-hit"    | cut -f 7 -d "," > blacklist/select_contaminant/contigs_to_keep.txt
tail -n +2 additional_data/blobtoolkit_GCA_000818735.3.csv | grep -e "Streptophyta" -e "no-hit" -v | cut -f 7 -d "," > blacklist/select_contaminant/contigs_to_remove.txt

# grep for contigs in vcf files
for l in $(seq 70 10 120); do 
   echo Checking vcf for read length ${l}
   grep -f blacklist/select_contaminant/contigs_to_remove.txt populations/populations_select_${l}_all_snps_whitelist/populations.snps.vcf    
done

# no contaminants found in assembled data

