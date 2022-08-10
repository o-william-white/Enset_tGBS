
library(phytools)
library(dplyr)

### read input trees

# input trees
t1 <- read.tree("best_trees/best.tre.fbp.raxml.support.rooted.newick")
t2 <- read.tree("best_trees/best.tre.iqtree.support.rooted.newick")

### change tree tip labels

# read sample metadata
sample.metadata <- read.csv("tGBS_metadata_phylogenetic_analysis.csv", colClasses = "character")

# sort sample meta data according to tip orders for t1 and t2
sample.metadata.t1 <- sample.metadata [ match(t1$tip.label, sample.metadata$popmap) , ]
sample.metadata.t2 <- sample.metadata [ match(t2$tip.label, sample.metadata$popmap) , ]

# check identical
print("Checking if tip labels and sample metadata are in the same order")
identical(t1$tip.label, sample.metadata.t1$popmap)
identical(t2$tip.label, sample.metadata.t2$popmap)

# change tip labels to reflect sample id and landrace name
t1$tip.label <- with(sample.metadata.t1, paste(sample_id, landrace, sep="_"))
t2$tip.label <- with(sample.metadata.t2, paste(sample_id, landrace, sep="_"))

### compare trees

# co-phylo object
obj <- cophylo(t1, t2)

### output

png(file="best_trees/compare_topology_raxml_vs_iqtree.png", width = 16*3, height = 23.7*3 , units = "cm", res=300)
plot(obj, fsize=0.75)
dev.off()

