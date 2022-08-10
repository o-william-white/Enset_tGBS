#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(vcfR)
library(adegenet)
library(poppr)
library(hierfstat)
library(cowplot)
library(ggplot2)
library(PerformanceAnalytics)
source("additional_scripts/HDplot.R")
source("additional_scripts/plot_hist_ggplot.R")

# input arguments
path_vcf = args[1]
path_md  = args[2]
out_prefix = args[3]

# eg 
# mkdir -p snp_summary_plots
# Rscript  additional_scripts/plot_snp_summary_stats.R  populations/populations_80_all_snps/populations.snps.vcf  tGBS_metadata_phylogenetic_analysis.csv  snp_summary_plots/summary

### read input data

# read vcf
vcf <- read.vcfR(path_vcf)

# read metadata
sample_metadata <- read.csv(path_md)

# basic checks
print("Check sample numbers per type")
table(sample_metadata$type)

print("Check vcf and metadata sample order is identical")
identical(as.character(colnames(vcf@gt)[-1]), 
          as.character(sample_metadata$sequence_id))


### per site stats

# convert vcf to genind
gi <- vcfR2genind(vcf, return.alleles=TRUE)

# add population info
gi@pop <- as.factor(sample_metadata$type)

# genind to hierfstat
hf <- genind2hierfstat(gi)

# calculate basic stats using hierfstat
bs <- basic.stats(hf)

# allele balance stats using hdplot
hd <- HDplot(vcf)

# extract depth data
dp <- extract.gt(vcf, element = "DP")

# calculate overall coverage (sum of depth) per site
co <- apply(dp, 1, function(x) sum(as.numeric(x), na.rm=T))

#  extract genotype data
gt <- extract.gt(vcf, element = "GT")

# calculate propotion of missing data per site
mi <- apply(gt, 1, function(x) sum(is.na(x))/length(x))

# dataframe of per site stats
per_site <- data.frame(
  ho = bs$perloc$Ho,
  fis = bs$perloc$Fis,
  ab = hd$ratio,
  co = co,
  mi = mi
)

# plot histograms
png(paste0(out_prefix, "_site_histograms.png"), height=8, width = 12, units="in", res=600)
plot_grid(plotlist = list(hist_ggplot(per_site$ho,  "Ho"),
                          hist_ggplot(per_site$fis, "Fis"),
                          hist_ggplot(per_site$ab,  "Allele balance"),
                          hist_ggplot(per_site$co,  "Coverage"),
                          hist_ggplot(per_site$mi,  "Missing", breaks = 10)),
          align = 'hv', axis = 'tblr', ncol = 3, nrow = 2)
dev.off()

# plot pairwise correlations
png(paste0(out_prefix, "_site_pairwise.png"), height=12, width = 12, units="in", res=600)
chart.Correlation(per_site)
dev.off()


### per indiv stats

# convert vcf to genind
gl <- vcfR2genlight(vcf)

# add population info
gl@pop <- as.factor(sample_metadata$type)

# bitwise genetic distance
bt <- bitwise.dist(gl, scale_missing = T)

# extract genotypes
gt <- extract.gt(vcf)

# proportion of heterozygous sites per indiv
# ignoring missing sites
mi <- apply(gt, 2, function(x) sum(x=="0/1", na.rm = T) / (length(x)-sum(is.na(x))))

# plot histograms
png(paste0(out_prefix, "_indiv_histograms.png"), height=4, width = 8, units="in", res=600)
plot_grid(plotlist = list(hist_ggplot(bt, breaks = 30, "Distance"),
                          hist_ggplot(mi, breaks = 30, "Proportion heterozygous sites")),
          align = 'hv', axis = 'tblr', ncol = 2)
dev.off()


