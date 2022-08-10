
# run as
# Rscript identify_high_depth_loci.R vcf output

# use the function commandArgs to pass commands to the script
args = commandArgs(trailingOnly=TRUE)

# get args
vcf <- args[1]
output <- args[2]

library(vcfR)

# read vcf
vcf <- read.vcfR(file = vcf)

# extract site depth
dp <- extract.gt(vcf, element = "DP", as.numeric = T)

# get site depth for each sample as a list 
list_snp_depth <- lapply(colnames(dp), function(x) dp[,x])
names(list_snp_depth) <- colnames(dp)

# function to plot histograms for each sample in list_snp_depth
plot_depth_hist <- function(x) {
  
  # test
  # x=1
  
  # get dat
  sample <- names(list_snp_depth)[x]
  snp_depth  <- list_snp_depth[[x]]
  
  # SNPs with NA depth values treated as zero 
  snp_depth [ is.na(snp_depth) ] <- 0  
    
  # loci names
  loci <- sapply(strsplit(names(snp_depth), ":"), getElement, 1)
  
  # loci names as.factor
  loci <- factor(loci, levels = unique(loci))
  
  # get max SNP depth per locus
  loci_max_depth <- tapply(snp_depth, loci, max)
  
  # get 95% quantile
  q <- quantile(loci_max_depth, probs = 0.95)
  
  # mean
  m <- mean(loci_max_depth)
  
  # n above 95% quantile
  n <- sum(loci_max_depth > q)
  
  hist(loci_max_depth, breaks=30, xlab = "Depth", main=sample)
  abline(v=m, col="blue", lty=1)
  abline(v=q, col="red", lty=2)
  legend("topright", inset = 0.025, legend = c(paste0("mean (", round(m, digits = 2), ")"), 
                                               paste0("q95 (", q, ")"),  
                                               paste(n, "loci with depth > q95")), 
         col=c("blue","red"), lty = c(1,2, 0))
  
}

# run function for first 9 samples and create a pdf
pdf(paste0(output, "/loci_depth_sample.pdf"), height = 12, width = 12)
par(mfrow=c(3,3))
lapply(1:9, plot_depth_hist)
dev.off()


# function to return loci with depth greater than 95% quantile of that sample
high_depth_loci <- function(snp_depth) {
  
  # test
  # snp_depth <- list_snp_depth[[1]]
 
  # SNPs with NA depth values treated as zero
  snp_depth [ is.na(snp_depth) ] <- 0

  # loci names
  loci <- sapply(strsplit(names(snp_depth), ":"), getElement, 1)
  
  # loci names as.factor
  loci <- factor(loci, levels = unique(loci))
  
  # get max SNP depth per locus
  loci_max_depth <- tapply(snp_depth, loci, max)
  
  # get 95% quantile
  q <- quantile(loci_max_depth, probs = 0.95)
  
  # loci with coverage more than 95% quantile
  high_covarage_loci <- names(loci_max_depth [ loci_max_depth > q  ])
  
  # return high_covarage_loci as factor
  factor(high_covarage_loci, levels = unique(loci))
  
}

# run function to find high depth loci for each sample
list_high_depth_loci <- lapply(list_snp_depth, high_depth_loci)

# unlist and count the number of samples in which each locus is identified as having a high depth
high_depth_loci_count <- table(unlist(list_high_depth_loci))

# identify loci that are identified as high coveriage in n sameples

# define threshold
threshold <- 25

# identify loci that are identified as high coveriage in 25 or more samples
loci_to_blacklist <- names(high_depth_loci_count) [ high_depth_loci_count > threshold ] 

# historgram showing the number of samples in which each locus is identified as having a high depth
# hist(high_depth_loci_count, breaks = 30, xlab="No. of samples in which each locus is identified as having a high depth", main="")
# abline(v=threshold, col="red", lty=2)

# write loci to blacklist
write.table(loci_to_blacklist,  paste0(output,"/blacklist.txt"), row.names = F, quote = F, col.names = F)

