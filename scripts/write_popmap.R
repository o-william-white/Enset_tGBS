
# read in metadata
sample.metadata <- read.csv("/data/scratch/mpx469/tGBS_enset_project/tGBS_metadata_phylogenetic_analysis.csv", header=TRUE, colClasses="character")

# create unique pop ids
population <- paste0("pop", 1:nrow(sample.metadata))

# cbind samples and population to create popmap
popmap <- cbind(sample.metadata$sequence_id, population)

# write popmap
write.table(file="popmap.txt", popmap, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

