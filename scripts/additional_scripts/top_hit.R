
# extract top hits from blastn output format 6
# requires dplyr

# run as
# Rscript top_hit.R <blastn_output>

# use the function commandArgs to pass commands to the script
args = commandArgs(trailingOnly=TRUE)

# load packages
library(dplyr)

# read input data
df <- read.table(args[1], header=FALSE)

# add column names 
colnames(df) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# select top hit, or first where the evalues are equal
th <- df %>% 
  group_by(qseqid) %>%
  slice_min(., order_by=evalue, with_ties=FALSE)

# write th
write.table(th, paste0(as.character(args[1]), "_top_hits"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

