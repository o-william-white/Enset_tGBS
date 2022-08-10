library(ggplot2)
library(cowplot)

dat <- read.table("gstacks/summary_gstacks.txt", header = T)

p1 <- ggplot(dat, aes(x=length, y=loci)) + 
  geom_bar(stat="identity") +
  scale_x_continuous(labels = seq(70,120,10), breaks = seq(70,120,10)) + 
  xlab("Sequence length") + 
  ylab("Loci") +
  theme_bw()

p2 <- ggplot(dat, aes(x=length, y=reads)) + 
  geom_bar(stat="identity") +
  scale_x_continuous(labels = seq(70,120,10), breaks = seq(70,120,10)) + 
  xlab("Sequence length") + 
  ylab("Reads") +
  theme_bw()

p3 <- ggplot(dat, aes(x=length, y=coverage)) + 
  geom_bar(stat="identity") +
  scale_x_continuous(labels = seq(70,120,10), breaks = seq(70,120,10)) + 
  xlab("Sequence length") + 
  ylab("Coverage") +
  theme_bw()

png("gstacks/summary_gstacks.png", width = 8, height = 8, res=300, units = "in")
plot_grid(plotlist = list(p1,p2,p3), align = "hv", axis = "tblr")
dev.off()
