
library(vcfR)
library(poppr)
library(hierfstat)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)

### read vcf

# read vcf file
vcf <- read.vcfR(paste0("populations/populations_select_80_all_snps_whitelist_blacklist_", args[1], "/populations.snps.vcf"))

### read sample metadata

# read data
sample_metadata <- read.csv("tGBS_metadata_population_genetics.csv", colClasses = "character")

# check numbers per type
# table(sample_metadata$type)

# check identical
print("Checking vcf and sample metadata names are identical")
identical(colnames(vcf@gt)[-1], sample_metadata$sequence_id)


### create genind object with population info

# convert vcf to genind
gi <- vcfR2genind(vcf, return.alleles=TRUE)

# check identical
print("Checking sample metadata and genind names are identical")
identical(sample_metadata$sequence_id, indNames(gi))

# add population info
gi@pop <- as.factor(sample_metadata$type)

# check numbers per type
# table(gi$pop)


### maf per population

# seppop
gi_sep <- seppop(gi, drop = FALSE)

# maf of each
maf <- lapply(gi_sep, function(x) as.numeric(minorAllele(x)))

# plot hist
par(mfrow=c(1,2))
hist(maf$Domesticated, main="Domesticated")
hist(maf$Wild,         main="Wild")

# maf equal to zero
#sum(maf$Domesticated == 0, na.rm = T)
#sum(maf$Wild == 0,         na.rm = T)

# wilcox test
# wilcox.test(maf$Domesticated, maf$Wild)

# mean values
# lapply(maf, function(x) mean(x, na.rm = T))


# plot histogram with separaate bar for zero values

# histogram breaks
hist_breaks <- seq(0.0, 0.5, by=0.05)

# add break to include maf of zero
hist_breaks <- c((hist_breaks[1] - 0.05), hist_breaks)

# check
# hist_breaks

# plot again - note the first bar equals for zero values
#par(mfrow=c(1,2))
#hist(maf$Domesticated, main="Domesticated", breaks = hist_breaks)
#hist(maf$Wild,         main="Wild", breaks = hist_breaks)

# get histogram data
hist_domesticated <- hist(maf$Domesticated, breaks = hist_breaks, plot = F) 
hist_wild         <- hist(maf$Wild,         breaks = hist_breaks, plot = F)

#hist_domesticated$counts
#hist_wild$counts

# function to fmt histogram data for ggplot
fmt_hist_data <-  function(x) {
  with(x, data.frame(mid     = mids,
                     low     = breaks[1:length(breaks)-1],
                     high    = breaks[2:length(breaks)],
                     density = density))
}

# frun function
hist_domesticated_df <- fmt_hist_data(hist_domesticated)
hist_wild_df         <- fmt_hist_data(hist_wild)

max_density <- 1.1*max(c(hist_domesticated_df$density),
                       c(hist_wild_df$density))

pal <- c("#0072B2", "#D55E00", "#E69F00", "black")
names(pal) <- c("Domesticated", "Feral", "Wild", "Outgroup")
#pal

#hist_domesticated_df
#hist_wild_df


# plot
plot_maf_density_domesticated <- plot_grid(
  
  ggplot(data = hist_domesticated_df[-1,], aes(x=mid, y=density)) + 
    geom_bar(stat = "identity", fill=pal["Domesticated"]) +
    scale_y_continuous(limits = c(0,max_density), expand = c(0,0)) +
    coord_flip() + 
    ggtitle("Domesticated") +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black", margin = margin(t=0,r=4,b=0,l=0, unit = "mm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5)),
  
  ggplot(data = hist_domesticated_df[1,], aes(x=mid, y=density)) + 
    geom_bar(stat = "identity", fill=pal["Domesticated"]) +
    scale_x_continuous(labels = c("0.0"), breaks = c(0), expand = c(0.5,0)) +
    scale_y_continuous(limits = c(0,max_density), expand = c(0,0)) +
    coord_flip() + 
    theme(axis.title = element_blank(),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black", margin = margin(t=0,r=4,b=0,l=0, unit = "mm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5)),
  
  ncol = 1, align = "hv", rel_heights = c(0.8,0.2))



plot_maf_density_wild <- plot_grid(
  
  ggplot(data = hist_wild_df[-1,], aes(x=mid, y=density)) + 
    geom_bar(stat = "identity", fill=pal["Wild"]) +
    coord_flip() +
    ggtitle("Wild") +
    scale_y_reverse(limits = c(max_density, 0), expand=c(0,0)) + 
    scale_x_continuous(position = "top") +
    theme(axis.title = element_blank(), 
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5)),
  
  ggplot(data = hist_wild_df[1,], aes(x=mid, y=density)) + 
    geom_bar(stat = "identity", fill=pal["Wild"]) +
    scale_x_continuous(breaks = c(0), position = "top", expand = c(0.5,0)) +
    scale_y_continuous(limits = c(0,max_density), expand = c(0,0)) +
    coord_flip() +
    scale_y_reverse(limits = c(max_density, 0), expand=c(0,0)) + 
    theme(axis.title = element_blank(), 
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5)),
  
  ncol = 1, align = "hv", rel_heights = c(0.8,0.2))


plot_overall <- plot_grid(plot_maf_density_wild,
                          plot_maf_density_domesticated,
                          rel_widths = c(0.48,0.52),     # to account for labels
                          align = "hv") 

# create output dir
dir.create(paste0("population_genetics/", args[1]), recursive = T, showWarnings = F)

# plot
png(paste0("population_genetics/", args[1], "/maf.png"), height = 3, width = 6, res = 600, units = "in")
plot_overall
dev.off()

# write data
write.table(hist_domesticated_df, paste0("population_genetics/", args[1], "/maf_domesticated.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
write.table(hist_wild_df,         paste0("population_genetics/", args[1], "/maf_wild.txt"),         col.names = T, row.names = F, quote = F, sep = "\t")

# check against base r
# hist(as.numeric(maf$Wild),         breaks = hist_breaks)
# hist(as.numeric(maf$Domesticated), breaks = hist_breaks) 


