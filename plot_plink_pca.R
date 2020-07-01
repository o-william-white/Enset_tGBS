# run as
# Rscript plot_plink_pca.R <plink_path>

# use the function commandArgs to pass commands to the script
args = commandArgs(trailingOnly=TRUE)

library(dplyr)

# read in pca data
eval <- read.table(file=paste0(args[1], "/plink.eigenval"), header=FALSE)
evec <- read.table(file=paste0(args[1], "/plink.eigenvec"), header=FALSE)

# set column names
colnames(evec) <- c("pop", "sample", paste0("e", 1:20))
colnames(eval) <- "eigenvalues"

# read in meta data with sample info
sample.metadata <- read.csv("/data/scratch/mpx469/tGBS_enset_project/tGBS_metadata.csv", header=TRUE, na.strings=NULL)

# remove data for "Disease" and "NA" samples
sample.metadata <- sample.metadata[sample.metadata$type != "Disease" & sample.metadata$type != "NA",]

# change levels of type
sample.metadata$type <- factor(sample.metadata$type, levels=c("Domestic", "Wild", "Outgroup"))

# remove endings from sample names in the pca data so they match meta data
evec$sample <- gsub(".unique.sorted", "", evec$sample)

# join dataframes to include all info
df <- left_join(evec, sample.metadata, by=c("sample" = "sequence_id"))

# filter to select only ingroup
df.in <- df %>% filter(type != "Outgroup")

# eigenvalues as percentages
# https://speciationgenomics.github.io/pca/
eval$pct <- eval$eigenvalues/sum(eval$eigenvalues)*100

pdf(paste0(args[1], "/plot_pca.pdf"))

plot(df$e1, df$e2, col=df$type, 
                         xlab=paste0("PC1 (", round(eval$pct[1],digits=2), " %)"), 
                         ylab=paste0("PC2 (", round(eval$pct[2],digits=2), " %)"))
legend("topright", inset=0.025, levels(df$type), col=1:3, pch=21)

dev.off()

