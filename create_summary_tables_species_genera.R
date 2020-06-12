
# get all files for species, genera
files.species <- list.files(path="summary_table_output/species/",     pattern = "species_count_*")
files.genera  <- list.files(path="summary_table_output/genera/",      pattern = "genera_count_*")

# read in all species and genera files as a list
list.files.species <- lapply(files.species, function(x) read.table(paste0("summary_table_output/species/", x),     header=TRUE, sep="\t"))
list.files.genera  <- lapply(files.genera,  function(x) read.table(paste0("summary_table_output/genera/", x),      header=TRUE, sep="\t"))

# set names, assuming files read in same order as specified
names(list.files.species) <- files.species
names(list.files.genera)  <- files.genera

# change column names for genera and species
for(i in files.species) {
  colnames(list.files.species[[i]]) <- c("species" ,i)
}

for(i in files.genera) {
  colnames(list.files.genera[[i]]) <- c("genera" ,i)
}

# join
df.species <- Reduce(function(...) merge(..., by='species', all.x=TRUE, all.y=TRUE), list.files.species)
df.genera  <- Reduce(function(...) merge(..., by='genera',  all.x=TRUE, all.y=TRUE), list.files.genera)

# change column names
colnames(df.species) <- gsub("species_count_", "", colnames(df.species))
colnames(df.genera)  <- gsub("genera_count_", "",  colnames(df.genera))

# write summary table
write.table(df.species, "summary_table_output/summary_species.txt",      col.names=TRUE, row.name=FALSE, quote=FALSE, sep="\t")
write.table(df.genera,  "summary_table_output/summary_genera.txt",       col.names=TRUE, row.name=FALSE, quote=FALSE, sep="\t")

q(save="no")

