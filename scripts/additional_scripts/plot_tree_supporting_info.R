
# load packages
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(ggnewscale)
library(cowplot)
library(tidyr)


### read trees

tre_raxml <- read.tree("best_trees/best.tre.tbe.raxml.support.rooted.newick")
tre_iqtre <- read.tree("best_trees/best.tre.iqtree.support.rooted.newick")

# plot
# plot(tre_raxml, cex=0.5)
# plot(tre_iqtre, cex=0.5)


### create pal for sample type

pal <- c("#0072B2", "#D55E00", "#E69F00", "black")
names(pal) <- c("Domesticated", "Semi-domesticated", "Wild", "Outgroup")
# pal

sha <- c(16,15,17,18)
names(sha) <- c("Domesticated", "Semi-domesticated", "Wild", "Outgroup")
# sha


### read and fmt sample metadata

# read sample metadata
sample_metadata <- read.csv("tGBS_metadata_phylogenetic_analysis.csv", colClasses = "character")

# change type to factor with levels
sample_metadata$type <- factor(sample_metadata$type, levels=c("Domesticated", "Semi-domesticated", "Wild", "Outgroup"))

# some minor changes to landrace names
sample_metadata$landrace <- gsub("Unnamed \\(semi-domesticated)", "Unnamed", sample_metadata$landrace)
sample_metadata$landrace <- gsub("horticultural", "hort.", sample_metadata$landrace)


### function to plot tree

plot_tree <- function(tre, format) {
  
  # test
  #tre = tre_raxml
  #format = "raxml"
  
  # tre = tre_iqtre
  # format = "iqtree"
  
  ### create ggtree object  
  
  # plot
  p <- ggtree(tre) %<+% sample_metadata 
  
  # add slot for mlg
  p$data$mlg <- rep(NA, nrow(p$data))
  
  
  
  ### add bootstap data
  
  if(format == "raxml") {
    # get bootstrap data for raxml
    bs <- filter(p$data, isTip != TRUE, label!="") %>%
      rename("bootstrap" = "label") %>%
      select(node, bootstrap, x, y) 
  } 
  
  
  if(format == "iqtree") {
    # get bootstap data for iqtree
    bs <- filter(p$data, isTip != TRUE, label!="") %>%
      rename("bootstrap" = "label") %>%
      select(node, bootstrap) %>%
      separate(bootstrap, into=c("uf_boot", "sh_alrt"), sep="/", remove = FALSE)
  }
  
  
  
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
  
  
  
  ### basic plots 
  
  if(format == "raxml") {
    
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
    
  } 
  
  
  if(format == "iqtree") {
    
    plot_tree <- p %<+% bs + 
      geom_tiplab(aes(label=paste(sample_id, landrace, sep=" ")), align = TRUE, size=3, offset = 0.00005) +
      geom_text2(aes(subset=!isTip & node!=0 & as.numeric(uf_boot) >= 95, label=bootstrap, x=branch), size=3, nudge_y = 0.75) +
      geom_tippoint(aes(col=type, shape=type), size = 3) + 
      geom_nodepoint(aes(subset=mlg=="mlg", x=x-0.015*max(p$data$x)), shape=8, size = 3) +
      scale_colour_manual(values=pal, na.translate=FALSE) + scale_shape_manual(values = sha) +  
      theme(legend.position = c(0.15, 0.75)) + 
      guides(linetype = "none", 
             colour = guide_legend("Type"),
             shape = guide_legend("Type")) + 
      xlim(0, 1.1*max(p$data$x))
    
  }

  return(plot_tree)

}

# run function
plot_tree_raxml <- plot_tree(tre_raxml, "raxml")
plot_tree_iqtre <- plot_tree(tre_iqtre, "iqtree")
 
# check output
# plot_tree_raxml
# plot_tree_iqtre


# pngs
png("best_trees/raxml_tbe.png", width = 16*3, height = 23.7*3, units = "cm", res = 600)
plot_tree_raxml
dev.off()

png("best_trees/iqtree.png", width = 16*3, height = 23.7*3, units = "cm", res = 600)
plot_tree_iqtre
dev.off()

