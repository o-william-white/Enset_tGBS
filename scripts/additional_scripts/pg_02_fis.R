
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



### low frequency alleles

# hist of maf
# hist(minorAllele(gi))

# filter for loci with maf >= 0.1
gi <- gi[loc=minorAllele(gi) >= 0.1,]

# check
# sum(minorAllele(gi) < 0.1)



### per locus statistics using hierfstat 

# convert genind to hierfstat fmt
hf <- genind2hierfstat(gi)

# calculate basic stats
bs <- basic.stats(hf)

# get fis
fis <- data.frame(bs$Fis)

# check
# head(fis)

# set breaks
breaks_hist <- seq(-1, 1, by=0.1)
# breaks_hist

# plot hist
hist(fis$Domesticated, xlab="", breaks = breaks_hist, main="Domesticated")
hist(fis$Wild,         xlab="", breaks = breaks_hist, main="Wild")

# get histogram data
hist_domesticated <- hist(fis$Domesticated, plot = F)
hist_wild         <- hist(fis$Wild,         plot = F)

# function to fmt histogram data for ggplot
fmt_hist_data <-  function(x) {
  with(x, data.frame(mid     = mids,
                     low     = breaks[1:length(breaks)-1],
                     high    = breaks[2:length(breaks)],
                     density = density))
}

# fmt data
hist_domesticated_df <- fmt_hist_data(hist_domesticated)
hist_wild_df         <- fmt_hist_data(hist_wild)

# check
# hist_domesticated_df
# hist_wild_df

# find max density for plot
max_density <- 1.1*max(c(hist_domesticated_df$density),
                       c(hist_wild_df$density))

pal <- c("#0072B2", "#D55E00", "#E69F00", "black")
names(pal) <- c("Domesticated", "Semi-domesitcated", "Wild", "Outgroup")
# pal

plot_fis_density_domesticated <- ggplot(data = hist_domesticated_df, aes(x=mid, y=density)) + 
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
        plot.title = element_text(hjust = 0.5))

plot_fis_density_wild <- ggplot(data = hist_wild_df, aes(x=mid, y=density)) + 
  geom_bar(stat = "identity", fill=pal["Wild"]) +
  coord_flip() +
  scale_y_reverse(limits = c(max_density, 0), expand=c(0,0)) + 
  scale_x_continuous(position = "top") +
  ggtitle("Wild") + 
  theme(axis.title = element_blank(), 
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))


plot_overall <- plot_grid(plot_fis_density_wild,
                          plot_fis_density_domesticated,
                          rel_widths = c(0.48,0.52),     # to account for labels
                          align = "hv") 

# create output dir
dir.create(paste0("population_genetics/", args[1]), recursive = T, showWarnings = F)

# plot
png(paste0("population_genetics/", args[1], "/fis.png"), height = 3, width = 6, res = 600, units = "in")
plot_overall
dev.off()

# write data
write.table(hist_domesticated_df, paste0("population_genetics/", args[1], "/fis_domesticated.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
write.table(hist_wild_df,         paste0("population_genetics/", args[1], "/fis_wild.txt"),         col.names = T, row.names = F, quote = F, sep = "\t")

