# install.packages('plyr', lib="/data/home/mpx469/software/R/3.6.1/", repos = 'https://cloud.r-project.org')

.libPaths("/data/home/mpx469/software/R/3.6.1/")
library(plyr)

# read in sample list
sample.list <- read.table("/data/scratch/mpx469/sample-list.txt", header=FALSE)

# read in count data
d2b.raw     <- read.table("count-data2bio-raw.txt",     header=FALSE, col.names="d2b.raw")
d2b.trimmed <- read.table("count-data2bio-trimmed.txt", header=FALSE, col.names="d2b.trimmed")
trimmomatic <- read.table("count-trimmomatic.txt",      header=FALSE, col.names="trimmomatic")

# plot histograms of reads per sample with a bin size of 5e5
pdf("plot-read-count-histograms.pdf")
par(mfrow=c(2,2))
hist(d2b.raw$d2b.raw,         breaks=seq(0, round_any(max(d2b.raw$d2b.raw),         5e5, f = ceiling),  by=5e5), main="Data2bio raw data" , xlab="Number of reads per sample")
hist(d2b.trimmed$d2b.trimmed, breaks=seq(0, round_any(max(d2b.trimmed$d2b.trimmed), 5e5, f = ceiling),  by=5e5), main="Data2bio trimmed data" , xlab="Number of reads per sample")
hist(trimmomatic$trimmomatic, breaks=seq(0, round_any(max(trimmomatic$trimmomatic), 5e5, f = ceiling),  by=5e5), main="Trimmomatic output" , xlab="Number of reads per sample")
dev.off()

# plot bar charts of read counts for first five samples 
df <- data.frame(d2b.raw, d2b.trimmed, trimmomatic)

row.names(df) <- sample.list$V1

df <- t(df)

pdf("plot-read-count-samples-1-5-barplot.pdf")
barplot(df[,1:5], cex.names=0.7, beside=TRUE, legend=rownames(df))
dev.off()

q(save="no")

