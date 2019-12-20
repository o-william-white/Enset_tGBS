# module add R/3.6.1
# R

# install packages as follows
# install.packages('plyr', lib="/data/home/mpx469/software/R/3.6.1/", repos = 'https://cloud.r-project.org')

.libPaths("/data/home/mpx469/software/R/3.6.1/")
library(plyr)
library(dplyr)
library(reshape2)

# read in sample list
sample.list <- read.table("/data/scratch/mpx469/sample-list.txt", header=FALSE)

# read in meta data
sample.metadata <- read.csv("/data/scratch/mpx469/GBS_metadata.csv", header=TRUE, na.strings=NULL)

# read in count data
raw         <- read.table("count-data2bio-raw.txt",     header=FALSE, col.names="raw")
trimmomatic <- read.table("count-trimmomatic.txt",      header=FALSE, col.names="trimmomatic")
cutadapt    <- read.table("count-cutadapt.txt",         header=FALSE, col.names="cutadapt")

# combine count data
df <- data.frame(raw, trimmomatic, cutadapt)

# add sample names
df$sample <- sample.list$V1

# caluclate the number removed by trimmomatic and cutadapt
df$pct.rm.by.trimmomatic <- (df$raw - df$trimmomatic) / df$raw
df$pct.rm.by.cutadapt <- (df$trimmomatic - df$cutadapt) / df$trimmomatic

# join sample metadata
df <- left_join(df, sample.metadata, by = c("sample" = "SAMPLE_ID"))%>%
      select(sample, raw, trimmomatic, pct.rm.by.trimmomatic, cutadapt, pct.rm.by.cutadapt, TYPE)

# remove NA and disease rows
df <- filter(df, !TYPE %in% c("NA", "Disease"))

# write summary table
write.table(df, "summary-read-counts.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# convert data to long format
df.long <- melt(df)

# plot histograms of reads per sample with a bin size of 5e5
pdf("plot-read-count-histograms-boxplot.pdf", height=8, width=8)
par(mfrow=c(2,2))
hist(df$raw,         breaks=seq(0, round_any(max(df$raw),         5e5, f = ceiling),  by=5e5), main="Raw data",            xlab="Number of reads per sample")
hist(df$trimmomatic, breaks=seq(0, round_any(max(df$trimmomatic), 5e5, f = ceiling),  by=5e5), main="Trimmomatic output" , xlab="Number of reads per sample")
hist(df$cutadapt,    breaks=seq(0, round_any(max(df$cutadapt),    5e5, f = ceiling),  by=5e5), main="Cutadapt output",     xlab="Number of reads per sample")
boxplot(value ~ variable, df.long, xlab = "", ylab = "Number of reads", main="Boxplot comparison")
dev.off()

png("plot-read-count-histograms-boxplot.png", height=8, width=8, units="in", res=300)
par(mfrow=c(2,2))
hist(df$raw,         breaks=seq(0, round_any(max(df$raw),         5e5, f = ceiling),  by=5e5), main="Raw data",            xlab="Number of reads per sample")
hist(df$trimmomatic, breaks=seq(0, round_any(max(df$trimmomatic), 5e5, f = ceiling),  by=5e5), main="Trimmomatic output" , xlab="Number of reads per sample")
hist(df$cutadapt,    breaks=seq(0, round_any(max(df$cutadapt),    5e5, f = ceiling),  by=5e5), main="Cutadapt output",     xlab="Number of reads per sample")
boxplot(value ~ variable, df.long, xlab = "", ylab = "Number of reads", main="Boxplot comparison")
dev.off()

