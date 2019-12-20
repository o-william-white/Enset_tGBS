
sample.metadata <- read.csv("/data/scratch/mpx469/GBS_metadata.csv", header=TRUE, na.strings=NULL)

# filter out disease and NA samples
sample.metadata <- sample.metadata[-which(sample.metadata$TYPE == "Disease" | sample.metadata$TYPE == "NA"),]

samples <- paste(sample.metadata$SAMPLE_ID, ".mapped.unique.sorted", sep="")

population.together <- rep(1, length(samples))
population.separate <- 1:length(samples)

popmap.together <- cbind(samples, population.together)
popmap.separate <- cbind(samples, population.separate)

write.table(file="popmap-selection-together.txt", popmap.together, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(file="popmap-selection-separate.txt", popmap.separate, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

q(save="no")

