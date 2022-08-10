
library(dplyr)

# read mlg data
mlg_farthest_popmap_corrected_single_rep <- read.table("estimate_mlg/mlg_farthest_bitwise_monophyletic_single_rep.txt", sep = "\t", header = TRUE, colClasses="character")

# write popmap for clone corrected dataset

# all samples
popmap_select_all <- mlg_farthest_popmap_corrected_single_rep %>%
  filter(type != "Semi-domesticated" & type != "Outgroup") %>%
  select(sequence_id, popmap)

# wild samples
popmap_select_wild <- mlg_farthest_popmap_corrected_single_rep %>%
  filter(type != "Semi-domesticated" & type != "Outgroup") %>%
  filter(type == "Wild") %>%
  select(sequence_id, popmap)

# domestic samples
popmap_select_domestic <- mlg_farthest_popmap_corrected_single_rep %>%
  filter(type != "Semi-domesticated" & type != "Outgroup") %>%
  filter(type == "Domesticated") %>%
  select(sequence_id, popmap)

print(paste("All = ", nrow(popmap_select_all)))
print(paste("Domesticated = ", nrow(popmap_select_domestic)))
print(paste("Wild = ", nrow(popmap_select_wild)))

# write popmaps
write.table(popmap_select_all,      "gstacks/popmap_select_all.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(popmap_select_wild,     "gstacks/popmap_select_wil.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(popmap_select_domestic, "gstacks/popmap_select_dom.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

# write sample metadata for population genetic dataset

# read sample metadata
sample_metadata <- read.csv("tGBS_metadata_phylogenetic_analysis.csv", colClasses="character")

# match popmap_all_select with sample metadata
sample_metadata_select <- sample_metadata[match(popmap_select_all$popmap, sample_metadata$popmap),]

# check identical
print("Checking sample metadata and popmap are in the same order")
identical(sample_metadata_select$popmap, popmap_select_all$popmap)

# write metadata
write.csv(sample_metadata_select, "tGBS_metadata_population_genetics.csv", quote = FALSE, row.names = FALSE)

