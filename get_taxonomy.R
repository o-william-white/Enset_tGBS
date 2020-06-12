# run as
# Rscript get_taxonomy.R <blastn_output>

# use the function commandArgs to pass commands to the script
args = commandArgs(trailingOnly=TRUE)

library(dplyr)



# read in ncbi rankedlineage.dmp
ncbi.lineage <- read.delim(file = "new_taxdump/rankedlineage_manual_edit.dmp",  header = FALSE, comment.char = "", quote = "", stringsAsFactors = FALSE, sep = "|")

# set colnames
colnames(ncbi.lineage) <- c("tax_id","tax_name","species","genus","family","order","class","phylum","kingdom","superkingdom")



# read blastn output
blast <- read.table(paste0("blastn_output/", args[1], "_blastn_out"), header=FALSE, sep = "\t")

# set column names for blast output
colnames(blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue","bitscore", "staxid", "ssciname")




# join taxonomy to blast output
blast <-  left_join(blast, ncbi.lineage, by=c("staxid"="tax_id")) 

# genus names incomplete for some taxa so create new genus column from tax_name
# filter(blast, genus_repaired == "")

# https://stackoverflow.com/questions/33683862/first-entry-from-string-split
blast$genus_repaired <- sapply(strsplit(blast$tax_name, split = " "), getElement, 1)

# some genera not to group due to additional symbols
# 'Nostoc
# [Enterobacter]
blast$genus_repaired <- gsub("'", "", blast$genus_repaired)
blast$genus_repaired <- gsub("\\[", "", blast$genus_repaired)
blast$genus_repaired <- gsub("\\]", "", blast$genus_repaired)

# get top sum_bitscore for each qseqid per taxid with filtered
taxid <- group_by(blast, qseqid, tax_name) %>%                          # group_by qseqid and staxids
            summarise(sum_bitscore = sum(bitscore)) %>%                 # sum bitscore per qseqid and staxids
            slice_max(n=1, order_by=sum_bitscore, with_ties=TRUE) %>%   # get top sum_bitscore for each qseqid (most likely taxonomy) with ties kept
            filter(as.logical(!anyDuplicated(qseqid)))                  # rm ties
 
# as above but for genera
genus <- group_by(blast, qseqid, genus_repaired) %>%
            summarise(sum_bitscore = sum(bitscore)) %>%
            slice_max(n=1, order_by=sum_bitscore, with_ties=TRUE) %>%
            filter(as.logical(!anyDuplicated(qseqid)))
 
# count per species and genus
count.taxid <- group_by(taxid, tax_name) %>% count(tax_name)
count.genus <- group_by(genus, genus_repaired) %>% count(genus_repaired)

# write species and genus counts
write.table(count.taxid, paste0("summary_table_output/species/species_count_", args[1]), col.names=TRUE, quote=FALSE, sep="\t", row.names=FALSE)
write.table(count.genus, paste0("summary_table_output/genera/genera_count_", args[1]), col.names=TRUE, quote=FALSE, sep="\t", row.names=FALSE)

