
library(vcfR)
library(poppr)
library(hierfstat)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
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


### exctract genotypes for domesticated and wild samples

# extract genotypes
geno <- extract.gt(vcf, element = 'GT')

# check
# geno[1:10,1:10]

# samples for domesticated and wild
samples_domesticated <- sample_metadata$sequence_id [ sample_metadata$type == "Domesticated" ]
samples_wild         <- sample_metadata$sequence_id [ sample_metadata$type == "Wild" ]

# subset geno for domesticated and wild
geno_domesticated <- geno [ , match(samples_domesticated, colnames(geno)) ]
geno_wild         <- geno [ , match(samples_wild,         colnames(geno)) ]


### calculate the proportion of heterozygous sites per individual

# function for proportion of heterozygous sites per individual
prop_het <- function(input_genotypes) {
  
  # count genotypes
  count_00 <- sum(input_genotypes == "0/0", na.rm = TRUE)
  count_01 <- sum(input_genotypes == "0/1", na.rm = TRUE)
  count_11 <- sum(input_genotypes == "1/1", na.rm = TRUE)
  
  # get sample size
  N <- count_00  + count_01 + count_11
  
  # retrn the proportion of heterozygous sites
  count_01 / N
  
}

# apply function to domesticated and wild smaples 
prop_het_domesticated <- apply(geno_domesticated, 2, prop_het)
prop_het_wild         <- apply(geno_wild,         2, prop_het)

# plot hists
# hist(prop_het_domesticated)
# hist(prop_het_wild)

# set breaks
breaks_hist <- breaks_hist <- seq(0, signif(max(c(prop_het_domesticated, prop_het_wild)),2)+0.01, by=0.01)
#breaks_hist

# set hist with set breaks
# hist(prop_het_domesticated, breaks = breaks_hist)
# hist(prop_het_wild,         breaks = breaks_hist)

# get histogram data
hist_domesticated <- hist(prop_het_domesticated, breaks = breaks_hist, plot = F)
hist_wild         <- hist(prop_het_wild,         breaks = breaks_hist, plot = F)

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

# ggplots
plot_het_density_domesticated <- ggplot(data = hist_domesticated_df, aes(x=mid, y=density)) + 
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

plot_het_density_wild <- ggplot(data = hist_wild_df, aes(x=mid, y=density)) + 
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

plot_overall <- plot_grid(plot_het_density_wild,
          plot_het_density_domesticated,
          rel_widths = c(0.48,0.52),     # to account for labels
          align = "hv") 

# create output dir
dir.create(paste0("population_genetics/", args[1]), recursive = T, showWarnings = F)

# plot
png(paste0("population_genetics/", args[1], "/prop_het.png"), height = 3, width = 6, res = 600, units = "in")
plot_overall
dev.off()

# write data
write.table(hist_domesticated_df, paste0("population_genetics/", args[1], "/prop_het_domesticated.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
write.table(hist_wild_df,         paste0("population_genetics/", args[1], "/prop_het_wild.txt"),         col.names = T, row.names = F, quote = F, sep = "\t")
