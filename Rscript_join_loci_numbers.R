# use the function commandArgs to pass commands to the script
args = commandArgs(trailingOnly=TRUE)

# set path to lib
.libPaths("/data/home/mpx469/software/R/3.6.1/")

# load packages
library(dplyr)

# get current dir
d <- getwd() 

# cd dir to input path
setwd(args[1])

# read in data

# duplicate sites
dupl <- read.table("duplicates", header=FALSE, col.names=c("count", "chr", "bp"))

# info
info <- read.table("loci_info", header=FALSE, col.names=c("locus", "chr", "bp"))

# join data
out <- left_join(dupl, info, by = c("chr"="chr", "bp"="bp"))

# write
write.table(out, "duplicates_loci", sep ="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

# return to base dir
setwd(d)

