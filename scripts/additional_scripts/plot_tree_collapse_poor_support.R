
# load packages
library(dplyr)
library(ggtree)
library(ape)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(ggnewscale)
library(cowplot)


### read tree

# read tre
tre <- read.tree("best_trees/best.tre.tbe.raxml.support.rooted.newick")

# collapse poorly supported nodes
tre <- as.polytomy(tre, feature = "node.label", fun=function(x) as.numeric(x) < 0.75)

# plot
# plot(tre, cex=0.5)



### create pal for sample type

pal <- c("#0072B2", "#D55E00", "#E69F00", "black")
names(pal) <- c("Domesticated", "Semi-domesticated", "Wild", "Outgroup")
#pal

sha <- c(16,15,17,18)
names(sha) <- c("Domesticated", "Semi-domesticated", "Wild", "Outgroup")
#sha



### read and fmt sample metadata

# read sample metadata
sample_metadata <- read.csv("tGBS_metadata_phylogenetic_analysis.csv", colClasses = "character")

# change type to factor with levels
sample_metadata$type <- factor(sample_metadata$type, levels=c("Domesticated", "Semi-domesticated", "Wild", "Outgroup"))

# some minor changes to landrace names
sample_metadata$landrace <- gsub("Unnamed \\(semi-domesticated)", "Unnamed", sample_metadata$landrace)
sample_metadata$landrace <- gsub("horticultural", "hort.", sample_metadata$landrace)



### create ggtree object  

# create ggtree object with group and sample metadata attached
p <- ggtree(tre) %<+% sample_metadata 

# add slot for mlg
p$data$mlg <- rep(NA, nrow(p$data))

# plot
# p + 
#   geom_tippoint(aes(col=type, shape=type)) + 
#   scale_colour_manual(values=pal, na.translate=FALSE) + 
#   scale_shape_manual(values = sha) + 
#   theme(legend.position = c(0.2,0.5)) +
#   guides(linetype = "none", 
#          colour = guide_legend("Type"), 
#          shape = guide_legend("Type"))



### add bootstap data

# filter for bs values more tahn or equal to 0.75
bs <- filter(p$data, isTip != TRUE, label!="") %>%
  rename("bootstrap" = "label") %>%
  select(node, bootstrap, x, y) 

# plot with bs data
# p %<+% bs + 
#   geom_text2(aes(subset=!isTip & node!=node.basal & as.numeric(bootstrap) >= 0.75, label=sprintf("%.2f", as.numeric(bootstrap)), x=branch), size=3, nudge_y = 0.75) +
#   geom_tippoint(aes(col=type, shape=type)) + 
#   scale_colour_manual(values=pal, na.translate=FALSE) + 
#   scale_shape_manual(values = sha) + 
#   theme(legend.position = c(0.2,0.5)) + 
#   guides(linetype = "none", 
#          colour = guide_legend("Type"), 
#          shape = guide_legend("Type"))


### read and fmt mlg

# read mlgs
mlg_farthest <- read.table("estimate_mlg/mlg_farthest_bitwise_monophyletic.txt", header = TRUE, col.names = c("sequence_id", "mlg_farthest", "popmap", "type"), colClasses = "character")

# check
# head(mlg_farthest)

# remove mlgs with only a single sample to reduce mlgs to colour
mlg_farthest <- group_by(mlg_farthest, mlg_farthest) %>% filter(n()!=1)

# as data.frame
mlg_farthest <- as.data.frame(mlg_farthest)

# get uniq mlgs after removing single sample mlgs
mlg_farthest_uniq <- unique(mlg_farthest$mlg_farthest)



### get mlgs as a list

# function to get popmap ids for a given mlg
myFunc <- function(mlg, popmap) {
  popmap$pop [ popmap$mlg == mlg  &  !is.na(popmap$mlg) ]
}

# for example
# myFunc(2, mlg_farthest)

# apply function to uniq mlgs with at least two sample
mlg_list <- lapply(mlg_farthest_uniq, function(x) myFunc(x, mlg_farthest))

# set names of list
names(mlg_list) <- paste0("mlg", mlg_farthest_uniq)

# all monophyletic
# table(unlist(lapply(mlg_list, function(x) is.monophyletic(tre, x))))

# mlg_list



### use list to find mrca for each mlg and to data

# change NA to "mlg" for mlg nodes  
for(i in 1:length(mlg_list)) {
  # get mlg node
  node_mlg <- ggtree::MRCA(p, mlg_list[[i]])
  # collapsed node set as TRUE for mlg_logical
  p$data$mlg [ p$data$node == node_mlg ] <- "mlg"
}



# add symbol to mlg branches
# p %<+% bs + 
#  geom_text2(aes(subset=!isTip & node!=0 & as.numeric(bootstrap) >= 0.75, label=sprintf("%.2f", as.numeric(bootstrap)), x=branch), size=3, nudge_y = 0.75) +
#  geom_tippoint(aes(col=type, shape=type)) + 
#  geom_nodepoint(aes(subset=mlg=="mlg", x=x-0.01*max(p$data$x)), shape=4, size=3) + 
#  scale_color_manual(values=pal, na.translate=FALSE) +
#  scale_shape_manual(values = sha) + 
#  theme(legend.position = c(0.15, 0.75)) + 
#  guides(colour=guide_legend("Type"),
#         shape=guide_legend("Type"),
#         linetype=FALSE)

### basic plots 
plot_tree <- p %<+% bs + 
  geom_tiplab(aes(label=paste(sample_id, landrace, sep=" ")), align = TRUE, size=3, offset = 0.00005) +
  geom_text2(aes(subset=!isTip & node!=0 & as.numeric(bootstrap) >= 0.75, label=sprintf("%.2f", as.numeric(bootstrap)), x=branch), size=3, nudge_y = 0.75) +
  geom_tippoint(aes(col=type, shape=type), size = 3) + 
  geom_nodepoint(aes(subset=mlg=="mlg", x=x-0.015*max(p$data$x)), shape=8, size = 3) +
  scale_colour_manual(values=pal, na.translate=FALSE) + scale_shape_manual(values = sha) +  
  theme(legend.position = c(0.15, 0.75)) + 
  guides(linetype = "none", 
         colour = guide_legend("Type"),
         shape = guide_legend("Type")) + 
  xlim(0, 1.1*max(p$data$x))

# plot_tree

png("best_trees/raxml_tbe_collapse_poor_supprt.png", width = 16*3, height = 23.7*3, units = "cm", res = 600)
plot_tree
dev.off()

