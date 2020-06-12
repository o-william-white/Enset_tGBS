
# get all files for species, genera
files.xcm <- list.files(path="summary_table_output/xcm_cov/", pattern = "*_cov")

# read in all species and genera files as a list
list.files.xcm <- lapply(files.xcm, function(x) read.table(paste0("summary_table_output/xcm_cov/", x), header=TRUE, sep="\t"))

# set names, assuming files read in same order as specified
names(list.files.xcm) <- gsub("_cov","", files.xcm)

# change column names
for(i in names(list.files.xcm)) {
  colnames(list.files.xcm[[i]]) <- c("chr" ,i)
}

# join dataframes in list
df.xcm <- Reduce(function(...) merge(..., by='chr', all.x=TRUE, all.y=TRUE), list.files.xcm)

# read in contig length
length.xcm <- read.table("contigs_xcm_lengths", col.names=c("chr","length"))

# sort length xcm to match df.xcm
length.xcm <- length.xcm[ match(df.xcm[,1], length.xcm[,1]) , ]

# get proportion of contigs recovered in blast search
df.xcm.prop <- cbind(df.xcm[,1], df.xcm[,-1] / length.xcm$length)

# set column name
colnames(df.xcm.prop)[1] <- "chr"

# write summary tables
write.table(df.xcm,      "summary_table_output/summary_xcm_bp.txt",            col.names=TRUE, row.name=FALSE, quote=FALSE, sep="\t")
write.table(df.xcm.prop, "summary_table_output/summary_xcm_bp_proportion.txt", col.names=TRUE, row.name=FALSE, quote=FALSE, sep="\t")

