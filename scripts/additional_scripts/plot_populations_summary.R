
library(ggplot2)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)

input_path = args[1]
output_path = gsub(".txt", ".png", input_path)

dat <- read.table(input_path, header = T)

p1 <- ggplot(dat, aes(x=length, y=loci)) + 
  geom_bar(stat="identity") +
  scale_x_continuous(labels = seq(70,120,10), breaks = seq(70,120,10)) + 
  xlab("Sequence length") + 
  ylab("Loci") +
  theme_bw()

p2 <- ggplot(dat, aes(x=length, y=sites)) + 
  geom_bar(stat="identity") +
  scale_x_continuous(labels = seq(70,120,10), breaks = seq(70,120,10)) + 
  xlab("Sequence length") + 
  ylab("Sites") +
  theme_bw()

p3 <- ggplot(dat, aes(x=length, y=filtered)) + 
  geom_bar(stat="identity") +
  scale_x_continuous(labels = seq(70,120,10), breaks = seq(70,120,10)) + 
  xlab("Sequence length") + 
  ylab("Filtered") +
  theme_bw()

p4 <- ggplot(dat, aes(x=length, y=variant)) + 
  geom_bar(stat="identity") +
  scale_x_continuous(labels = seq(70,120,10), breaks = seq(70,120,10)) + 
  xlab("Sequence length") + 
  ylab("Variant") +
  theme_bw()

png(output_path, width = 8, height = 8, res=300, units = "in")
plot_grid(plotlist = list(p1,p2,p3,p4), align = "hv", axis = "tblr")
dev.off()
