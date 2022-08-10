
library(dplyr)
library(tidyr)
library(ggplot2)
library(ape)
library(ggtree)
library(aplot)
library(reshape2)
library(RColorBrewer)

### 4) EDIT THE PATH TO YOUR COPY of FinestructureLibrary.R
source("additional_scripts/FinestructureLibrary.R", chdir = TRUE) # read in the R functions, which also calls the needed packages

# read coancestry output
chunkfile <- read.table("radpainter/populations.haps.clone.correct_chunks.out", header = TRUE, row.names = 1)

# see histogram of coancestry
# hist(melt(chunkfile)$value, main="", xlab="Coancestry")

# cap max conacestry as 110
# chunkfile [ chunkfile > 110 ] <- 110

# rowname to column
chunkfile$sequnce_id <- row.names(chunkfile)
row.names(chunkfile) <- NULL

# melt to long
chunkfile <- melt(chunkfile)

# change column names
colnames(chunkfile)[1:2] <- c("x", "y")

chunkfile$value [ chunkfile$x == chunkfile$y ] <- NA

# check
# head(chunkfile)

# read coancestry tree
treexml<-xmlTreeParse("radpainter/populations.haps.clone.correct_chunks.mcmcTree.xml") ## read the tree as xml format
ttree<-extractTree(treexml) ## extract the tree into ape's phylo format

# plot(ttree)

# read mlg
mlg <- read.table("estimate_mlg/mlg_farthest_bitwise_monophyletic_single_rep.txt", header = TRUE)

# as factor
mlg$type <- factor(mlg$type, levels = c("Domesticated", "Semi-domesticated", "Wild", "Outgroup"))

# create pal for sample type
pal <- c("#0072B2", "#D55E00", "#E69F00")
names(pal) <- c("Domesticated", "Semi-domesticated", "Wild")
#pal

# shape
sha <- c(16,15,17)
names(pal) <- c("Domesticated", "Semi-domesticated", "Wild")

# get ggtree object

# for ploting with tip labels
# p <- ggtree(tre, ladderize = FALSE) + geom_tiplab(size=3, align = TRUE) 
# p <- p + xlim(0, 1.2*max(p$data$x))

# without tip labels
p <- ggtree(ttree, ladderize = FALSE) %<+% mlg 

l <- max(p$data$x)

p <- p + geom_tiplab(aes(label = ""), align = TRUE) +
   geom_tippoint(aes(colour = type, shape = type, x = l), size = 1) + 
   scale_colour_manual(values = pal) +
   guides(colour = guide_legend(title = "Type", 
                                override.aes = list(size = 4),
                                ncol = 1,
                                direction = "vertical",
                                label.position = "top",
                                title.position = "left",
                                label.theme = element_text(vjust = 0.5, size = 12, angle=90)),
          shape = guide_legend("Type")) +
   theme(legend.title = element_text(size = 12, angle = 90), 
         legend.position = "right")
    
# p

#### plot heat map coancestry againsts the tree

# set levels to match tree tip labels
chunkfile$x <- factor(chunkfile$x, levels = ttree$tip.label)
chunkfile$y <- factor(chunkfile$y, levels = ttree$tip.label)


### 5) EXECUTE THE CODE ABOVE AND THE REST OF THE CODE BELOW
## make some colours
some.colors<-MakeColorYRP() # these are yellow-red-purple
some.colorsEnd<-MakeColorYRP(final=c(0.2,0.2,0.2)) # as above, but with a dark grey final for capped values

heatmap_coancestry <- ggplot(chunkfile, aes(x, y, fill=value)) + 
   geom_tile() +
   scale_fill_gradientn(colours = some.colorsEnd, na.value = "transparent") + 
   labs(x = "", y = "", fill = "Coancestry") +
   theme_bw() + 
   theme(legend.position = "right",
         legend.title = element_text(size = 12, angle = 90),
         axis.text.x = element_text(size = 4, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
         axis.text.y = element_text(size = 4, colour = "black"),
         axis.ticks = element_blank(),
         panel.grid = element_blank(),
         legend.text = element_text(size = 8)) + 
   guides(fill = guide_colourbar(#direction = "horizontal", 
                                 title.position = "left",
                                 barwidth = 1.5, 
                                 barheight = 8))

# heatmap_coancestry

heatmap_coancestry_no_labels <- ggplot(chunkfile, aes(x, y, fill=value)) + 
   geom_tile() +
   scale_fill_gradientn(colours = some.colorsEnd, na.value = "transparent") + 
   labs(x = "", y = "", fill = "Coancestry") +
   theme_bw() + 
   theme(legend.position = "right",
         legend.title = element_text(size = 12, angle = 90),
         axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks = element_blank(),
         panel.grid = element_blank(),
         legend.text = element_text(size = 8)) + 
   guides(fill = guide_colourbar(#direction = "horizontal", 
      title.position = "left",
      barwidth = 1.5, 
      barheight = 8))

# heatmap_coancestry_no_labels

# heatmap_coancestry %>% insert_left(p, width = 0.2) 


# plot with tree
# dev.off()
pdf("radpainter/plot_coancestry_radpainter.pdf", width = 10, height = 8)
heatmap_coancestry_no_labels %>% insert_left(p, width = 0.2) 
dev.off()

