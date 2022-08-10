
### load libraries
library(vcfR)
library(poppr)
library(ggplot2)
library(reshape2)


### read vcf

# read vcf file
vcf <- read.vcfR("populations/populations_80_all_snps_clean/populations.snps.vcf")


### read and format sample metadata

# read sample metadata
sample_metadata <- read.csv("tGBS_metadata_phylogenetic_analysis.csv", colClasses = "character")

# change type to factor with levels
sample_metadata$type <- factor(sample_metadata$type, levels=c("Domesticated", "Semi-domesticated", "Wild", "Outgroup"))


### farthest mlg threshold

farthest <- read.table("estimate_mlg/farthest_threshold.txt")$V1
# farthest


### get genlight 

# convert vcf to genlight
gl <- vcfR2genlight(vcf)

# sample metadata and gl in same order
print("Checking sample metadata and genlight object are in the same order")
identical(indNames(gl), sample_metadata$sequence_id)

# add population info
gl@pop <- sample_metadata$type


### get snpclone

# convert to snpclone
sc <- as.snpclone(gl)


### function to return subsamples of population <pop>, of size <n>, repeated <rep> times 

subsample_pop <- function(sc, pop, n, rep) {
  # sample names for pop
  sample_names <- sc$ind.names [ sc$pop == pop ]
  # subsample pop for n samples over rep
  subsample <- replicate(rep, sample(sample_names, n, replace = F), simplify = F)
  # return list
  return(subsample)
}

# subsample_test <- subsample_pop(sc=sc, pop="Domesticated", n=5, rep=10)
# subsample_test


### function to count the number of mlgs in a subsample of a given population based on a specified distance threshold

count_mlg <- function(sc, subsample, thres) {
  # filter sc for subsample
  sc <- sc [ sc$ind.names %in% subsample , ]
  # mlg filter based on threshold
  mlg.filter(sc, bitwise.dist, scale_missing = T) <- thres
  # mlg for pop
  mlg <- mlg.table(sc, plot = F)
  # return number of domesticated mlgs
  return(length(mlg))
}

# count_mlg(sc=sc, subsample=subsample_test[[1]], thres=farthest)
# sapply(subsample_test, function(x) count_mlg(sc=sc, subsample=x, thres=farthest))

# set sample sizes
# sum(sc$pop == "Domesticated")
sample_sizes <- seq(from=10, to=220, by=10)
print("Using the following sample sizes")
sample_sizes

# run, this takes some time...
print("Starting rarefaction")
list_mlg <- lapply(sample_sizes, function(n) {
  subsample_n <- subsample_pop(sc=sc, pop="Domesticated", n=n, rep=100)
  sapply(subsample_n, function(x) count_mlg(sc, x, farthest))
})

# as matrix
matrix_mlg <- do.call("cbind", list_mlg)
colnames(matrix_mlg) <- sample_sizes
rownames(matrix_mlg) <- paste0("repeat_", 1:nrow(matrix_mlg))
# matrix_mlg

# write csv
write.csv(matrix_mlg, "estimate_mlg/mlg_accumulation.csv", quote = F, row.names = F)

# melt to long format
mlg_long <- melt(matrix_mlg, varnames = c("Repeat", "Sample_size"), value.name = "MLGs")

# write csv
write.csv(mlg_long, "estimate_mlg/mlg_accumulation_long.csv", quote = F, row.names = F)

