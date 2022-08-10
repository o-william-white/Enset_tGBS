
# load packages
library(dplyr)
library(ape)
library(phylosignal)
library(phylobase)
library(PerformanceAnalytics)


### read and fmt tree

# read tre
tre <- read.tree("best_trees/best.tre.tbe.raxml.support.rooted.newick")

# remove node bs labels or else phylobase throws loads of errors when we run phylo4
tre$node.label <- NULL

# plot
# plot(tre, cex=0.5)


### read sample metadata

# read data
sample_metadata <- read.csv("tGBS_metadata_phylogenetic_analysis.csv", colClasses = "character")


### drop tips to select the domesticated samples

# sort metadata to match tree
sample_metadata <- sample_metadata [ match(tre$tip.label, sample_metadata$popmap) , ]

# check identical
print("Checking sample metadata and tree tip labels are in the same order")
identical(sample_metadata$popmap, tre$tip.label)

# drop tips for those not domestic
tre_domestic <- drop.tip(tre, tre$tip.label [ sample_metadata$type != "Domesticated" ])

# check
# plot(tre_domestic)


### read trait data

# read data
trait_data <- read.csv("additional_data/trait_datasetzonal_filter_corrected.csv")

# check
# head(trait_data)

# popmap to rownames
rownames(trait_data) <- trait_data$popmap
trait_data <- trait_data[,-1]

# drop observation
trait_data <- trait_data[,-1]

# check 
# head(trait_data)


### test for phylogentic signal using phylosignal package

# create phylo4d object
ps <- phylo4d(x=tre_domestic, tip.data=trait_data)

# add positive and negative control
tipData(ps)$random <- rnorm(length(tre_domestic$tip.label)) 
tipData(ps)$BM     <- rTraitCont(as(ps, "phylo"))

# measures of signal
ps_out <- phyloSignal(ps)
# ps_out

# write result
write.table(ps_out$stat,   "phylosignal/phylosig_stats.txt",  sep = "\t", col.names=TRUE, quote=FALSE)
write.table(ps_out$pvalue, "phylosignal/phylosig_pvalue.txt", sep = "\t", col.names=TRUE, quote=FALSE)

# plot traits
#barplot.phylo4d(ps,     tree.type = "phylo", tree.ladderize = FALSE, center = FALSE, scale = FALSE, bar.lwd = 0.2, tip.font = 1)
#dotplot.phylo4d(ps,     tree.type = "phylo", tree.ladderize = FALSE, center = FALSE, scale = FALSE, bar.lwd = 0.2, tip.font = 1)
#gridplot.phylo4d(ps,    tree.type = "phylo", tree.ladderize = FALSE, center = FALSE, scale = FALSE, bar.lwd = 0.2, tip.font = 1)


### locate signal with lipa

# lipa
#lipa_ps <- lipaMoran(ps)
#lipa_ps$p.adjust <- apply(lipa_ps$p.value, 2, function(x) p.adjust(x, "BH"))

# here, we use a proximity matrix based on the number of nodes to ignore the effect of long terminal branches and focus on clades
#lipa_ps     <- lipaMoran(p4d, prox.phylo = "nNodes")

# phylosig plot
#barplot.phylo4d(ps,     bar.col=(lipa_ps$p.value < 0.05) + 1,     center = TRUE , scale = TRUE, bar.lwd = 5)
#barplot.phylo4d(ps,     bar.col=(lipa_ps$p.adjust < 0.05) + 1,     center = TRUE , scale = TRUE, bar.lwd = 5)

# write pvalues
# write.table(lipa_ps$p.value,     "results_table_lipa_pvalues.txt",  quote = FALSE, row.names = TRUE, sep="\t")

# write pdf of plot
#dev.off()
#pdf("phlosignal_lipa.pdf", width = 6.29921*2, height = 9.330709*3.7)
#barplot.phylo4d(ps,     bar.col=(lipa_ps$p.value < 0.05) + 1,     center = TRUE , scale = TRUE, bar.lwd = 5)
#barplot.phylo4d(ps,     bar.col=(lipa_ps$p.adjust < 0.05) + 1,     center = TRUE , scale = TRUE, bar.lwd = 5)
#dev.off()




