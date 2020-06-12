
# read in metadata
sample.metadata <- read.csv("/data/scratch/mpx469/tGBS_enset_project/tGBS_metadata.csv", header=TRUE, na.strings=NULL)

# filter out disease and NA samples
sample.metadata <- sample.metadata[-which(sample.metadata$type == "Disease" | sample.metadata$type == "NA"),]

# add samtools suffix to sequence id
samples <- paste(sample.metadata$sequence_id, ".unique.sorted", sep="")

# create unique pop ids
population <- paste0("pop", 1:length(samples))

# cbind samples and population to create popmap
popmap <- cbind(samples, population)

# write popmap
write.table(file="popmap.txt", popmap, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

