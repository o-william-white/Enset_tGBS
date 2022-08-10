
# load packages
library(ape)
library(dplyr)
library(ggplot2)
library(ggtree)


### read trees

tree_raxml <- read.tree("best_trees/best.tre.fbp.raxml.support.rooted.newick")
tree_iqtre <- read.tree("best_trees/best.tre.iqtree.support.rooted.newick")

# plot
# plot(tree_raxml, cex=0.5)
# plot(tree_iqtre, cex=0.5)


### read and fmt sample metadata

# read sample metadata
sample_metadata <- read.csv("tGBS_metadata_phylogenetic_analysis.csv", colClasses = "character")

# change type to factor with levels
sample_metadata$type <- factor(sample_metadata$type, levels=c("Domesticated", "Semi-domesticated", "Wild", "Outgroup"))


### read and fmt mlg

# read mlgs
mlg_farthest <- read.table("estimate_mlg/mlg_farthest_bitwise.txt", header = TRUE, col.names = c("sequence_id", "mlg_farthest"), colClasses = "character")

# join popmap info and sample type from metadata
mlg_farthest_popmap <- left_join(mlg_farthest, sample_metadata, by="sequence_id") %>%
  select(sequence_id, mlg_farthest, popmap, type)

# function to get popmap ids for a given mlg
myFunc <- function(mlg, popmap) {
  popmap$pop [ popmap$mlg == mlg  &  !is.na(popmap$mlg) ]
}

# for example
# myFunc(1, mlg_farthest_popmap)
# myFunc(2, mlg_farthest_popmap)

# lapply to get popmaps ids for each mlg as list
mlg_farthest_list <- lapply(unique(mlg_farthest_popmap$mlg_farthest), function(x) myFunc(x, mlg_farthest_popmap))

# set names
names(mlg_farthest_list) <- paste0("mlg_", unique(mlg_farthest_popmap$mlg_farthest))

# check
# mlg_farthest_list

# subset list for mlgs with more than one sample
mlg_farthest_list_multiple <- mlg_farthest_list [ lapply(mlg_farthest_list, length) > 1 ] 

# check again
# mlg_farthest_list_multiple

# how many mlgs
print(paste(length(mlg_farthest_list), "MLGs identified"))
# how many mlgs with more than one sample
print(paste(length(mlg_farthest_list_multiple), "MLGs have more than one sample"))


### test monophyly of each mlg with the raxml tree 

# not all monophyletic!
monophyly_df <- data.frame(raxml = unlist(lapply(mlg_farthest_list_multiple, function(x) is.monophyletic(tree_raxml, x))),
                           iqtre = unlist(lapply(mlg_farthest_list_multiple, function(x) is.monophyletic(tree_iqtre, x))))
           
# identify non-monophyletic mlgs
non_monophyly_df <- monophyly_df [ monophyly_df$raxml == FALSE | monophyly_df$iqtre == FALSE , ]

print(paste(nrow(non_monophyly_df), "MLGs are non-monophyletic"))
print(non_monophyly_df)

# get non-monophyletic mlgs names
non_monophyly_names <- rownames(non_monophyly_df)


### plot all mlgs individually


# create ggtree objects
p_raxml <- ggtree(tree_raxml, layout = "circular") %<+% sample_metadata
p_iqtre <- ggtree(tree_iqtre, layout = "circular") %<+% sample_metadata

# raxml

# create emply plot list
plot_list_raxml <- vector("list", length = length(mlg_farthest_list_multiple))

# copy
p_raxml_i <- p_raxml

# add slot for point col
p_raxml_i$data$point_col <- rep("x", nrow(p_raxml_i$data))

# for loop
for(i in 1:length(mlg_farthest_list_multiple)) {
  # set point col for mlg
  p_raxml_i$data$point_col [ match(mlg_farthest_list_multiple[[i]], p_raxml_i$data$label) ] <- "mlg"
  # plot
  p_raxml_i2 <- p_raxml_i + 
    geom_tippoint(aes(col=point_col), size=2) + 
    geom_tiplab(aes(label=landrace), align = TRUE) +
    theme(legend.position = "none") + 
    ggtitle(names(mlg_farthest_list_multiple)[i])
  # add plot to list  
  plot_list_raxml[[i]] = p_raxml_i2
  # reset point col
  p_raxml_i$data$point_col <- rep("x", nrow(p_raxml_i$data))
}

# plot raxml tree with each mlg individually
pdf("estimate_mlg/plot_trees_raxml_mlg_farthest.pdf", height = 14, width = 14)
for (i in 1:length(mlg_farthest_list_multiple)) {
  print(plot_list_raxml[[i]])
}
dev.off()


# iqtree

# create emply plot list
plot_list_iqtre <- vector("list", length = length(mlg_farthest_list_multiple))

# copy
p_iqtre_i <- p_iqtre

# add slot for point col
p_iqtre_i$data$point_col <- rep("x", nrow(p_iqtre_i$data))

# for loop
for(i in 1:length(mlg_farthest_list_multiple)) {
  # set point col for mlg
  p_iqtre_i$data$point_col [ match(mlg_farthest_list_multiple[[i]], p_iqtre_i$data$label) ] <- "mlg"
  # plot
  p_iqtre_i2 <- p_iqtre_i + 
    geom_tippoint(aes(col=point_col), size=2) + 
    geom_tiplab(aes(label=landrace), align = TRUE) +
    theme(legend.position = "none") +
    ggtitle(names(mlg_farthest_list_multiple)[i])
  # add plot to list  
  plot_list_iqtre[[i]] = p_iqtre_i2
  # reset point col
  p_iqtre_i$data$point_col <- rep("x", nrow(p_iqtre_i$data))
}


# plot iqtre tree with each mlg individually
pdf("estimate_mlg/plot_trees_iqtree_mlg_farthest.pdf", height = 14, width = 14)
for (i in 1:length(mlg_farthest_list_multiple)) {
  print(plot_list_iqtre[[i]])
}
dev.off()


### change mlgs for samples which are non monophyletic

# filter monophyletic
mlg_farthest_popmap_mono <-      filter(mlg_farthest_popmap, !mlg_farthest %in% gsub("mlg_", "", non_monophyly_names))

# filter non-monophyletic
mlg_farthest_popmap_non_mono <-  filter(mlg_farthest_popmap, mlg_farthest %in% gsub("mlg_", "", non_monophyly_names))

# see non mono
# mlg_farthest_popmap_non_mono

# change non mono mlg to separate mlgs
mlg_farthest_popmap_non_mono$mlg_farthest <- as.character(max(as.numeric(mlg_farthest_popmap_mono$mlg_farthest))+1:nrow(mlg_farthest_popmap_non_mono) )

# check 
# mlg_farthest_popmap_non_mono

# recombine mono and corrected  mlgs
mlg_farthest_popmap_corrected <- rbind(mlg_farthest_popmap_mono,
                                       mlg_farthest_popmap_non_mono)

# compare the number of mlgs before and after
# length(unique(mlg_farthest_popmap$mlg_farthest))
# length(unique(mlg_farthest_popmap_corrected$mlg_farthest))

# difference of 9
# nrow(mlg_farthest_popmap_non_mono) - length(non_monophyly_names) 


### get a single representative per mlg 

# select first duplicate of each mlg
mlg_farthest_popmap_corrected_single_rep <- mlg_farthest_popmap_corrected [ !duplicated(mlg_farthest_popmap_corrected$mlg_farthest) , ]

# check that the number of mlg retained matched the exp
# nrow(mlg_farthest_popmap_corrected_single_rep)
# length(unique(mlg_farthest_popmap_corrected$mlg_farthest))

# 250 samples
# nrow(mlg_farthest_popmap_corrected)

# 153 mlgs
# nrow(mlg_farthest_popmap_corrected_single_rep)
# length(unique(mlg_farthest_popmap_corrected$mlg_farthest))

# see mlg data
# head(mlg_farthest_popmap_corrected)
# head(mlg_farthest_popmap_corrected_single_rep)


# write mono mlgs
write.table(mlg_farthest_popmap_corrected, "estimate_mlg/mlg_farthest_bitwise_monophyletic.txt", col.names = c("sequence_id", "mlg_farthest", "popmap", "type"), row.names = FALSE, quote=FALSE, sep="\t")

# write mono mlgs with single representative per mlg
write.table(mlg_farthest_popmap_corrected_single_rep, "estimate_mlg/mlg_farthest_bitwise_monophyletic_single_rep.txt", col.names = c("sequence_id", "mlg_farthest", "popmap", "type"), row.names = FALSE, quote=FALSE, sep="\t")

# for corrected mlgs - lapply to get popmaps ids for each mlg as list
mlg_farthest_corrected_list <- lapply(unique(mlg_farthest_popmap_corrected$mlg_farthest), function(x) myFunc(x, mlg_farthest_popmap_corrected))

# histogram of mlg size
# length(mlg_farthest_corrected_list)

# how many mlgs with more than one sample
# sum(unlist(lapply(mlg_farthest_corrected_list, function(x) length(x) > 1)))

# histogram of mlg disitribution
# hist(unlist(lapply(mlg_farthest_corrected_list, length)), xlab="MLG size", main="")

# table(unlist(lapply(mlg_farthest_corrected_list, length)))

