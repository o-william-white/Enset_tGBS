
library(reshape2)

s <- read.table("summary_gstacks")

png("summary_gstacks.png", width = 7, height = 7, res=300, units = "in")

par(mfrow=c(2,2))

barplot(s[,"loci"], names.arg = row.names(s), col="gray80", ylab="Loci", cex.axis = 0.6, cex.names = 0.6, cex.lab=0.6)
title(main="A", cex.main=0.75, adj = 0)

barplot(s[,"reads"], names.arg = row.names(s), col="gray80", ylab="Reads", cex.axis = 0.6, cex.names = 0.6, cex.lab=0.6)
title(main="B", cex.main=0.75, adj = 0)

barplot(s[,"coverage"], names.arg = row.names(s), col="gray80", ylab="Coverage", cex.axis = 0.6, cex.names = 0.6, cex.lab=0.6)
title(main="C", cex.main=0.75, adj = 0)

dev.off()

