
library(dplyr)
library(stringr)

# read sample metadata
sample_metadata <- read.csv("tGBS_metadata_phylogenetic_analysis.csv", colClasses = "character")

# select sequence_id and type column
sample_metadata <- sample_metadata [ , c("sequence_id", "type") ]

# create "by_species"
sample_metadata$by_sample <- paste(sample_metadata$sequence_id, sample_metadata$type, sep = "_")

# change new label for outgroup samples to "Outgroup"
sample_metadata$by_sample [ str_detect(sample_metadata$by_sample, "Outgroup") ] <- "Outgroup"

# change "type" to "by_pop"
colnames(sample_metadata)[2] <- "by_pop"

# check
# head(sample_metadata)

# read mlg
mlg <- read.table("estimate_mlg/mlg_farthest_bitwise_monophyletic_single_rep.txt", header = TRUE)

# for samples not defined as single representatives of a mlg - change "by_pop" to "xxx"
sample_metadata$by_pop [ !sample_metadata$sequence_id %in% mlg$sequence_id ] <- "xxx"

# for samples not defined as single representatives of a mlg - change "by_sample" to "xxx"
sample_metadata$by_sample [ !sample_metadata$sequence_id %in% mlg$sequence_id ] <- "xxx"

# check 
# sample_metadata

# write sets
write.table(sample_metadata[,c(1,2)], "dsuite/by_pop.txt",    col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(sample_metadata[,c(1,3)], "dsuite/by_sample.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

