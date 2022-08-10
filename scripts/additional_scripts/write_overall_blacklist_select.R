# run as
# Rscript write_overall_blacklist.R length

# use the function commandArgs to pass commands to the script
args = commandArgs(trailingOnly=TRUE)

# get args
length <- args[1]

# read individual blacklists
blastn         <- read.table(paste0("blacklist/select_blastn/blastn_",                 length, "/blacklist.txt"),  col.names = "loci")
high_depth     <- read.table(paste0("blacklist/select_high_depth/high_depth_",         length, "/blacklist.txt"),  col.names = "loci")
duplicate_loci <- read.table(paste0("blacklist/select_duplicate_loci/duplicate_loci_", length, "/blacklist.txt"),  col.names = "loci")
  
# get unique loci for blacklist
blacklist_overall <- unique(c(blastn$loci, high_depth$loci, duplicate_loci$loci))
  
# summary dataframe
summary_df <- t(data.frame(
  blastn          = nrow(blastn),
  high_depth      = nrow(high_depth),
  duplicate_loci  = nrow(duplicate_loci),
  overall         = length(blacklist_overall)
))
  
# write blacklist
write.table(blacklist_overall, paste0("blacklist/select_overall_blasklist/blacklist_", length, ".txt"), row.names = F, col.names = F, quote = F)
  
# write blacklist
write.table(summary_df, paste0("blacklist/select_overall_blasklist/summary_", length, ".txt"),  row.names = T, col.names = F, quote = F, sep = "\t")

