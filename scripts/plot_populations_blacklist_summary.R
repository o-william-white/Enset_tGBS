
library(reshape2)

s <- read.table("summary_single_snp_blacklist")[,c(1,4)]
a <- read.table("summary_all_snps_blacklist")[,c(1,4)]

s <- t(s)
a <- t(a)

pdf("summary_populations_blacklist.pdf", width = 7, height = 4)

par(mfrow=c(1,2))

barplot(s, beside=TRUE, col=c("gray40","gray80"), main="single snp", cex.axis = 0.5, cex.names = 0.5, cex.main = 0.75, adj = 0)
legend("topright", inset=0.025, legend = c("loci", "snps"), fill=c("gray40","gray80"), cex = 0.5)

barplot(a, beside=TRUE, col=c("gray40","gray80"), main="all snps", cex.axis = 0.5, cex.names = 0.5, cex.main = 0.75, adj = 0)
legend("topleft", inset=0.025, legend = c("loci", "snps"), fill=c("gray40","gray80"), cex = 0.5)

dev.off()

