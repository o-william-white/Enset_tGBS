
library(vcfR)
library(poppr)
library(reshape2)

### set up dir

if (!dir.exists("estimate_mlg")){
  dir.create("estimate_mlg")
}

### read vcf

# read vcf file
vcf <- read.vcfR("populations/populations_80_all_snps_clean/populations.snps.vcf")


### read and format sample metadata

# read sample metadata
sample_metadata <- read.csv("tGBS_metadata_phylogenetic_analysis.csv", colClasses = "character")

# change type to factor with levels
sample_metadata$type <- factor(sample_metadata$type, levels=c("Domesticated", "Semi-domesticated", "Wild", "Outgroup"))


### estimate multilocus genotypes

# convert vcf to genlight
gl <- vcfR2genlight(vcf)

# sample metadata and gl sample order identical
print("Checking sample metadata and genlight object have the same sample order")
identical(as.character(gl@ind.names), as.character(sample_metadata$sequence_id))

# add population info
gl@pop <- sample_metadata$type

# convert to snpclone
sc <- as.snpclone(gl)

# histograms of pairwise genetic distances for all and each population
# bin sizes done by eye to be broadly similar
png("estimate_mlg/histograms_pairwise_genetic_distance.png", height = 8, width = 8, units = "in", res = 600)
par(mfrow=c(2,2))
hist(bitwise.dist(sc,                              scale_missing = TRUE), breaks=50, xlim=c(0,0.33), xlab="Distance", main = "All samples")
hist(bitwise.dist(popsub(sc, "Domesticated"),      scale_missing = TRUE), breaks=40, xlim=c(0,0.33), xlab="Distance", main = "Domesticated")
hist(bitwise.dist(popsub(sc, "Semi-domesticated"), scale_missing = TRUE), breaks=20, xlim=c(0,0.33), xlab="Distance", main = "Semi-domesticated")
hist(bitwise.dist(popsub(sc, "Wild"),              scale_missing = TRUE), breaks=30, xlim=c(0,0.33), xlab="Distance", main = "Wild")
dev.off()

# filter stats 
fs <- filter_stats(sc, bitwise.dist, scale_missing = TRUE, plot=TRUE)

# see https://cran.rstudio.com/web/packages/poppr/vignettes/mlg.html#custom-custom 
# section "Choosing a threshold"

# set threshold between the first small peak and second larger peak
# the initial peak likely represents clones differentiated by a small set of random mutations
# this isn't obvious for our data

# another method is to look for the largest gap between all putative thresholds
# using the cutoff_predictor() fuction 

# three possible algorithms
# nearest neighbor
# farthest neighbor (default)
# average neighbor (UPGMA)

# get thresholds
farthest <- cutoff_predictor(fs$farthest$THRESHOLDS)
average  <- cutoff_predictor(fs$average$THRESHOLDS)
nearest  <- cutoff_predictor(fs$nearest$THRESHOLDS)

# add to plot
abline(v=farthest, col="#E41A1C", lty=3)
abline(v=average,  col="#377EB8", lty=4)
abline(v=nearest,  col="#4DAF4A", lty=5)

png("estimate_mlg/filter_stats_bitwise.png")
filter_stats(sc, bitwise.dist, scale_missing = TRUE, plot=TRUE)
abline(v=farthest, col="#E41A1C", lty=3)
abline(v=average,  col="#377EB8", lty=4)
abline(v=nearest,  col="#4DAF4A", lty=5)
dev.off()

# create new sc object for each threshold
sc_farthest <- sc
sc_average  <- sc
sc_nearest  <- sc

# filter sc based on each threshold
mlg.filter(sc_farthest, bitwise.dist, scale_missing = TRUE) <- farthest 
mlg.filter(sc_average,  bitwise.dist, scale_missing = TRUE) <- average
mlg.filter(sc_nearest,  bitwise.dist, scale_missing = TRUE) <- nearest

# create counts, vectors, and matrices of multilocus genotypes
#mlg.table(sc_farthest)
#mlg.table(sc_average)
#mlg.table(sc_nearest)

# create data frame with samples and mlg
mlg_farthest <- data.frame(sample=sc_farthest@ind.names,
                           mlg=mlg.vector(sc_farthest))

mlg_average  <- data.frame(sample=sc_average@ind.names,
                           mlg=mlg.vector(sc_average))

mlg_nearest <- data.frame(sample=sc_nearest@ind.names,
                           mlg=mlg.vector(sc_nearest))

# order by mlg
mlg_farthest <- mlg_farthest[order(mlg_farthest$mlg),]
mlg_average  <- mlg_average[order(mlg_average$mlg),]
mlg_nearest  <- mlg_nearest[order(mlg_nearest$mlg),]

# write mlg
write.table(mlg_farthest, "estimate_mlg/mlg_farthest_bitwise.txt", quote=FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)
write.table(mlg_average,  "estimate_mlg/mlg_average_bitwise.txt",  quote=FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)
write.table(mlg_nearest,  "estimate_mlg/mlg_nearest_bitwise.txt",  quote=FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)

# pairwise genetic distances as long df
bitwise_dist <- as.matrix(bitwise.dist(gl, scale_missing = TRUE))
bitwise_dist_long <- melt(bitwise_dist, varnames = c("sampleA", "sampleB"), value.name = "bitwise_distance")

# pairwise genetic distances as long df for lower triangle
bitwise_dist_lower_tri <- bitwise_dist
bitwise_dist_lower_tri [ lower.tri(bitwise_dist_lower_tri) ] <- NA
bitwise_dist_lower_tri <- melt(bitwise_dist_lower_tri, varnames = c("sampleA", "sampleB"), value.name = "bitwise_distance")
bitwise_dist_lower_tri$bitwise_distance [ bitwise_dist_lower_tri$sampleA == bitwise_dist_lower_tri$sampleB ] <- NA
bitwise_dist_lower_tri <- bitwise_dist_lower_tri [ !is.na(bitwise_dist_lower_tri$bitwise_distance) , ]

# compare plots 
# hist(bitwise_dist_lower_tri$bitwise_distance, breaks = 50)
# hist(bitwise.dist(gl, scale_missing = T), breaks = 50)

# write table with pairwise distances 
write.table(bitwise_dist_long,      "estimate_mlg/pairwise_bitwise_distances_long.txt",           col.names = T, row.names = F, sep = "\t", quote = F)
write.table(bitwise_dist_lower_tri, "estimate_mlg/pairwise_bitwise_distances_long_lower_tri.txt", col.names = T, row.names = F, sep = "\t", quote = F)

# write threshold
write.table(farthest, "estimate_mlg/farthest_threshold.txt", col.names = F, row.names = F)

