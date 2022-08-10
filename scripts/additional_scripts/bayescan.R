
library(vcfR)
library(dplyr)
library(adegenet)
source("additional_scripts/plot_R.r")

#### check results

# read dat
sel1=read.table("bayescan/bayescan_output_1/populations.snps.maf.bayescan.sel",colClasses="numeric")
sel2=read.table("bayescan/bayescan_output_2/populations.snps.maf.bayescan.sel",colClasses="numeric")
sel3=read.table("bayescan/bayescan_output_3/populations.snps.maf.bayescan.sel",colClasses="numeric")

# plot
pdf("bayescan/convergence_logL.pdf")
par(mfrow=c(3,1))
with(sel1, plot(logL, type="l"))
with(sel2, plot(logL, type="l"))
with(sel3, plot(logL, type="l"))
dev.off()

# read dat
dat1 <- read.table("bayescan/bayescan_output_1/populations.snps.maf.bayescan_fst.txt", header = T)
dat2 <- read.table("bayescan/bayescan_output_2/populations.snps.maf.bayescan_fst.txt", header = T)
dat3 <- read.table("bayescan/bayescan_output_3/populations.snps.maf.bayescan_fst.txt", header = T)

# default plot
pdf("bayescan/outliers.pdf", height = 12, width = 12)
par(mfrow=c(2,2))
outliers1 <- plot_bayescan(dat1, FDR = 0.05)
outliers2 <- plot_bayescan(dat2, FDR = 0.05)
outliers3 <- plot_bayescan(dat3, FDR = 0.05)
dev.off()

# common elements of multiple lists
# example
# a <- c(1,3,5,7,9)
# b <- c(3,6,8,9,10)
# c <- c(2,3,4,5,7,9)
# Reduce(intersect, list(a,b,c))

# get index of outliers across all three runs
index_outliers <- Reduce(intersect, list(outliers1$outliers,
                                         outliers2$outliers,
                                         outliers3$outliers))


### get outlier info

# read vcf
vcf <- read.vcfR("bayescan/populations.snps.sorted.dw.maf.recode.vcf")

# get outlier chr, pos, id and CLocus
vcf_dat <- data.frame(Chrom     = getCHROM(vcf),
                      Pos       = getPOS(vcf),
                      Id        = getID(vcf),
                      Ref       = getREF(vcf),
                      Alt       = getALT(vcf),
                      Site      = paste(getID(vcf), getALT(vcf), sep = "."),
                      CLocus    = paste0("CLocus_", sapply(strsplit(getID(vcf),":"), getElement, 1)))

# filter vcf data for outlier loci
outlier_dat <- vcf_dat[index_outliers,]

# print number of outliers
print(paste(nrow(outlier_dat), "outlier snps across", length(unique(outlier_dat$CLocus)), "loci"))

# a positive value of alpha suggests diversifying selection, whereas negative values suggest balancing or purifying selection
# dat1$alpha [ index_outliers ]
# dat2$alpha [ index_outliers ]
# dat3$alpha [ index_outliers ]


### are outliers simply due to noise

# see distribution of fst
# hist(dat1$qval)

col <- rep(1, nrow(dat1))
col[ index_outliers ] <- 2

# q-value, which is the FDR analogue of the p-value

# a q-value, which is the FDR analogue of the p-value 
# (note that a q-value is only defined in the context of multiple testing, 
# whereas a p-value is defined on a single test). The q-value of given locus
# is the minimum FDR at which this locus may become significant.

#dat1$minus_log10_qval <- -log10(dat1$qval)
#plot(dat1$minus_log10_qval, ylab="-log10(qval)", col = col)

pdf("bayescan/outliers_qval_fst.pdf", height = 10, width = 10)
par(mfrow=c(2,1))
plot(dat1$qval, ylab="qval", col = col)
abline(h=0.05, col="blue", lty=2)
plot(dat1$fst, ylab="fst", col = col)
dev.off()


### create bed files 

# bed file for outliers 
outlier_bed <- data.frame(Chrom = outlier_dat$Chrom,
                          Start = outlier_dat$Pos-1,
                          End   = outlier_dat$Pos)

# bed file for all sites
vcf_bed <- data.frame(Chrom = vcf_dat$Chrom,
                      Start = vcf_dat$Pos-1,
                      End   = vcf_dat$Pos)

# write bed files
write.table(outlier_bed, "bayescan/outliers.bed",   col.names = F, row.names = F, quote = F, sep = "\t")
write.table(vcf_bed,     "bayescan/background.bed", col.names = F, row.names = F, quote = F, sep = "\t")



### get allele frequencies in the domesticated and wild populations

# create genind object
gi <- vcfR2genind(vcf, return.alleles = T)

# read population definitions used to format the bayescan data
pop_def <- read.table("bayescan/vcf_to_bayescan_definitions.txt", col.names = c("sequence_id", "pop"), colClasses = "character")

# check names are identical
print("Checking sample names are identical")
identical(indNames(gi), pop_def$sequence_id)

# set population id
gi$pop <- factor(pop_def$pop)

# genind as genpop
gp <- genind2genpop(gi)

# allele freq as data.frame
allele_freq <- data.frame(t(tab(gp)))

# add rownames to column
allele_freq <- cbind(Site = row.names(allele_freq), allele_freq)
row.names(allele_freq) <- NULL

# get ID as factor
allele_freq$Site <- as.character(allele_freq$Site)
allele_freq$Id <- factor(sapply(strsplit(allele_freq$Site, "\\."), "[", 1))

# check
# head(allele_freq)

# get allele counts per population in format ref/alt 
Domesticated_ref_alt <- tapply(allele_freq$Domesticated, allele_freq$Id, function(x) paste(x, collapse = ","))
Wild_ref_alt         <- tapply(allele_freq$Wild,         allele_freq$Id, function(x) paste(x, collapse = ","))

# check 
# head(Domesticated_ref_alt)
# head(Wild_ref_alt)



### add allele frequencies to outlier data

outlier_dat$Domesticated_ref_alt <- sapply(outlier_dat$Id, function(x) Domesticated_ref_alt[x])
outlier_dat$Wild_ref_alt         <- sapply(outlier_dat$Id, function(x) Wild_ref_alt[x])


# check
# head(outlier_dat)


### get snpeff annotation

# read annotation
#snp_ann <- read.table("../snpeff/snpsift_wide_ann_fmt.txt", header = T, sep = "\t")
#head(snp_ann)

### add annotations to outlier data

#outlier_dat <- left_join(outlier_dat, snp_ann, by="Site")


### write outlier data

# select relevant columns
outlier_dat2 <- outlier_dat %>% select(Chrom, Pos, Ref, Alt, CLocus, Domesticated_ref_alt, Wild_ref_alt)

# remove annotation columns that are all NA
# outlier_dat2 <- outlier_dat2 [ -which(apply(outlier_dat2, 2, function(x) all(is.na(x)))) ] 

# change NA to ""
# outlier_dat2 [ is.na(outlier_dat2) ] <- ""

# check again
# head(outlier_dat2)

# write table
write.table(outlier_dat2, "bayescan/outlier_info.txt", col.names = T, sep="\t", quote = F)


