
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
dupl <- read.table("out-duplicated-sites.txt", header=FALSE, col.names=c("count", "chr", "bp"))

# info
info <- read.table("out-sumstats-info.txt", header=FALSE, col.names=c("locus", "chr", "bp"))

# join data
out <- left_join(dupl, info, by = c("chr"="chr", "bp"="bp"))

# write
write.table(out, "out-duplicated-sites-and-sumstats-info.txt", sep ="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

# return to base dir
setwd(d)

