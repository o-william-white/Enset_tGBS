
library(ape)
library(dplyr)
library(stringr)

# tree
tre <- read.tree("best_trees/best.tre.tbe.raxml.support.rooted.newick")

# plot(tre)

# only want a single outgroup tip - drop two of three outgroups
tre <- drop.tip(tre, "pop82")
tre <- drop.tip(tre, "pop162")

# check again
# plot(tre)

# read sample metadata
sample_metadata <- read.csv("tGBS_metadata_phylogenetic_analysis.csv", colClasses = "character")

# create new label 
sample_metadata$new_label <- paste(sample_metadata$sequence_id, sample_metadata$type, sep = "_")

# get df with popmap and new label
df <- select(sample_metadata, popmap, new_label)

# change outgroup label
df$new_label [ str_detect(df$new_label, "Outgroup") ] <- "Outgroup"

# change tip labels
tre$tip.label <- df$new_label [ match(tre$tip.label, df$popmap) ]

# plot(tre)

# remove node labels
tre$node.label <- NULL

# write tree
write.tree(tre, "dsuite/raxml.renamed.newick")

