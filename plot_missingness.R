
# use the function commandArgs to pass commands to the script
args = commandArgs(trailingOnly=TRUE)

imiss <- read.table(paste0("missingness_", args[1], "_", args[2], "_output/missingness_", args[1], "_", args[2], ".imiss"), header=TRUE)
lmiss <- read.table(paste0("missingness_", args[1], "_", args[2], "_output/missingness_", args[1], "_", args[2], ".lmiss"), header=TRUE)

pdf(file=paste0("missingness_", args[1], "_", args[2], "_output/missingness_", args[1], "_", args[2], ".pdf"), height=4, width=8)
par(mfrow=c(1,2))
with(lmiss, hist(F_MISS, breaks=80, xlab="Proportion missing", main="Loci", cex.lab=0.6, cex.axis=0.6, cex.main=0.6))
with(imiss, hist(F_MISS, breaks=80, xlab="Proportion missing", main="Individuals", cex.lab=0.6, cex.axis=0.6, cex.main=0.6))
dev.off()

