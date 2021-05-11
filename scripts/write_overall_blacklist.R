
# function to read in and summarise blacklist info
summarise_blacklists <- function(length, dataset) {
  
  # read individual blacklists
  blastn         <- read.table(paste0("/data/scratch/mpx469/tGBS_enset_project/blacklists/blastn/blastn_", length, "_", dataset, "/blastn_blacklist.txt"),                         col.names = "loci")
  loci_depth     <- read.table(paste0("/data/scratch/mpx469/tGBS_enset_project/blacklists/loci_depth/loci_depth_", length, "_", dataset, "_output/loci_depth_blacklist.txt"),          col.names = "loci")
  duplicate_loci <- read.table(paste0("/data/scratch/mpx469/tGBS_enset_project/blacklists/duplicate_loci/duplicate_loci_", length, "_", dataset, "_output/duplicate_loci_blacklist.txt"),  col.names = "loci")
  
  # get unique loci for blacklist
  blacklist_overall <- unique(c(blastn$loci, loci_depth$loci, duplicate_loci$loci))
  
  # summary dataframe
  summary_df <- t(data.frame(
    blastn          = nrow(blastn),
    loci_depth      = nrow(loci_depth),
    duplicate_loci  = nrow(duplicate_loci),
    overall         = length(blacklist_overall)
  ))
  
  # write blacklist
  write.table(blacklist_overall, paste0("blacklist_", length, "_", dataset, ".txt"), row.names = F, col.names = F, quote = F)
  
  # write blacklist
  write.table(summary_df, paste0("summary_", length, "_", dataset, ".txt"),  row.names = T, col.names = F, quote = F, sep = "\t")
  
}

# run function for each dataset
for(length in seq(70,120, by=10)){
  for(dataset in c("single_snp", "all_snps")) {
    summarise_blacklists(length, dataset)
  }
}

