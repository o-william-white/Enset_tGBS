# read in data
total <- read.table("flagstat-plot/reads-total.txt", header=FALSE)
unique <- read.table("flagstat-plot/reads-unique.txt", header=FALSE)
unique.mapped <- read.table("flagstat-plot/reads-unique-mapped.txt", header=FALSE)

# rbind data in long format
df <- rbind(total, unique, unique.mapped)

# add class info
df <- cbind(df, c(rep("total", nrow(total)), rep("unique", nrow(unique)), rep("unique mapped", nrow(unique))))

colnames(df) <- c("count", "class")

pdf(file="flagstat-plot/plot-flagstat.pdf")
boxplot(count ~ class, df, xlab = "", sub=paste("Mean proportion of uniquely mapped reads = ", signif(mean(unique.mapped$V1 / unique$V1), digits=2), "%", sep =""))
dev.off()

png(file="flagstat-plot/plot-flagstat.png")
boxplot(count ~ class, df, xlab = "", sub=paste("Mean proportion of uniquely mapped reads = ", signif(mean(unique.mapped$V1 / unique$V1), digits=2), "%", sep =""))
dev.off()

q(save="no")


