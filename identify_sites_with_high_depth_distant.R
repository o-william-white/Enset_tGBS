# identify sites with consistently high coverage
# requires dplyr

# run as
# Rscript identify_sites_with_high_depth.R <LENGTH> <DATASET>

# use the function commandArgs to pass commands to the script
args = commandArgs(trailingOnly=TRUE)

library(dplyr)



# get path to dir with site depths output from vcftools
dir.path <- paste0("site_depth_distant_", args[1], "_", args[2], "_output/")



# read loci info (required later to join chrom pos to loci no)
loci.info <- read.table(paste0("loci_info_distant_", args[1], "_", args[2]), col.names = c("loci", "chr", "pos"))

# paste chr and pos in new column called chr.pos
loci.info$chr.pos <- paste(loci.info$chr, loci.info$pos)



# list files for each sample
files <- list.files(path = dir.path , pattern = ".ldepth")

# read files into list (only need columns 1:3)
list.depth <- lapply(files, function(x) read.table(paste0(dir.path, x), header=TRUE)[,1:3])

# set names
names(list.depth) <- gsub(".ldepth", "", files)

# take a sample of first nine samples
list.depth.sample <- list.depth[1:9]



# function to plot histograms for first nine samples 
plotSample <- function(x) {
  
  # mean
  m <- mean(list.depth[[x]]$SUM_DEPTH)
  
  # 95% quantile
  q <- quantile(list.depth[[x]]$SUM_DEPTH, probs = 0.95)
  
  # n above 95% quantile
  n <- nrow(list.depth[[x]][ list.depth[[x]]$SUM_DEPTH > q ,])
  
  hist(list.depth[[x]]$SUM_DEPTH, breaks=100, xlab = "Depth", main=x)
  abline(v=q, col="red", lty=2)
  abline(v=m, col="blue", lty=1)
  legend("topright", inset = 0.025, legend = c(paste0("mean (", round(m, digits = 2), ")"), paste0("q95 (", q, ")"),  paste(n, "> q95")), col=c("blue","red"), lty = c(1,2, 0))
  
}



# run function and plot first nine samples
pdf(paste0("site_depth_distant_" , args[1], "_", args[2], "_sample.pdf"), height = 12, width = 12)
par(mfrow=c(3,3))
lapply(names(list.depth.sample), plotSample)
dev.off()



# function to filter for sites depth > 95% quantile in each sample
filterSiteDepth <- function(x) {
  
  # 95% quantile
  q <- quantile(list.depth[[x]]$SUM_DEPTH, probs = 0.95)
  
  # select site with coverage more than 95% quantile
  list.depth[[x]] <- list.depth[[x]] [ list.depth[[x]]$SUM_DEPTH > q , ]
  
}


# run function
list.depth.filtered <- lapply(names(list.depth), filterSiteDepth)

# set names 
names(list.depth.filtered) <- gsub(".ldepth", "", files)


# function to get unique chrom pos combinations for high depth sites
uniqChromPos <- function(x) {
  
  paste(list.depth.filtered[[x]]$CHROM, list.depth.filtered[[x]]$POS)  
  
}


# run function
list.depth.filtered.uniq <- lapply(names(list.depth.filtered), uniqChromPos)

# unlist and count unique chrom pos for sites across all samples
site.counts <- table(unlist(list.depth.filtered.uniq))

# identify chr and pos of sites with high depth in > 3 samples (just over 10% of 28 samples) 
chr.pos <- names(site.counts[ site.counts > 3 ])

# change into dataframe
chr.pos <- data.frame(chr.pos)

# join df to loci.info from populations sumstats output to get loci numbers
df <- left_join(chr.pos, loci.info, by="chr.pos")

# get uniq loci for blacklist
blacklist <- select(df, loci) %>% unique

# write blacklist
write.table(blacklist, paste0("blacklist_site_depth_distant_", args[1], "_", args[2]), row.names = FALSE, col.names = FALSE, quote = FALSE)


