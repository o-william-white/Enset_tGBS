# read in input data
reads.raw         <- read.table(file="read_length.EXS11ID000851.digested.txt", header=FALSE)
reads.trim        <- read.table(file="read_length.EXS11ID000851.digested.trimmed.txt", header=FALSE)
reads.trimmomatic <- read.table(file="read_length.EXS11ID000851.digested.trimmomatic.txt", header=FALSE)

# create line plot
pdf(file="plot.read_length.EXS11ID000851.digested.pdf")
plot(reads.raw$V2,reads.raw$V1,type="l",xlab="read length",ylab="occurences",col="blue", main="EXS11ID000851")
lines(reads.trim$V2,reads.trim$V1, col="red")
lines(reads.trimmomatic$V2,reads.trimmomatic$V1, col="green")
legend(min(reads.raw$V2), max(reads.raw$V1), legend=c("raw", "trimmed", "trimmomatic"), col=c("blue", "red", "green"), lty=c(1,1,1), ncol=1)
dev.off()

png(file="plot.read_length.EXS11ID000851.digested.png")
plot(reads.raw$V2,reads.raw$V1,type="l",xlab="read length",ylab="occurences",col="blue", main="EXS11ID000851")
lines(reads.trim$V2,reads.trim$V1, col="red")
lines(reads.trimmomatic$V2,reads.trimmomatic$V1, col="green")
legend(min(reads.raw$V2), max(reads.raw$V1), legend=c("raw", "trimmed", "trimmomatic"), col=c("blue", "red", "green"), lty=c(1,1,1), ncol=1)
dev.off()

q(save="no")
