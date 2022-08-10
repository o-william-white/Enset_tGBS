
library(dplyr)

# read sample metadata
md <- read.csv("tGBS_metadata.csv")

# drop SEQUENCE_ID
md <- md[,-which(colnames(md) == "SEQUENCE_ID")]

# change SAMPLE_ID to sequence_id
colnames(md)[which(colnames(md) == "SAMPLE_ID")] <- "sequence_id"

# drop POPULATION
md <- md[,-which(colnames(md) == "POPULATION")]

# change TYPE to type
colnames(md)[which(colnames(md) == "TYPE")] <- "type"

# replace "Domestic" with "Domesticated"
md$type <- gsub("Domestic", "Domesticated", md$type, fixed = T)

# replace "Domesticated_collection" with "Domesticated"
md$type <- gsub("Domesticated_collection", "Domesticated", md$type)

# vector for semi-domesticated samples
samples_semi_domesticated <- c("P3S62ID554", "P3S61ID553", "P2S77ID344", "P2S76ID345", "EXS17ID000914","EXS16ID000913", "P3S60ID912")

# change type for semi-domesticated samples
md$type [ match(samples_semi_domesticated, md$sequence_id) ] <- "Semi-domesticated"

# P1EN082 Eventricosum_Maurelii
# orginally identified as outgroup
# md$type [ md$sequence_id == "P1EN082" ]
# This sample groups with domesticates
# changed TYPE to Domesticated
md$type [ md$sequence_id == "P1EN082" ] <- "Domesticated"

# P1EN079 Eventricosum
# orginally identified as outgroup
# md$type [ md$sequence_id == "P1EN079" ]
# This sample groups with domesticates, but is of unknown provenance
# changed type to NA 
md$type [ md$sequence_id == "P1EN079" ] <- NA

# P1EN080 Eventricosum
# P1EN081 Eventricosum
# These samples group with other species of Ensete (sister to E. glaucum), but are of unknown provenance
# md$type [ md$sequence_id == "P1EN080" ]
# md$type [ md$sequence_id == "P1EN081" ]
# changed type to NA
md$type [ md$sequence_id == "P1EN080" ] <- NA
md$type [ md$sequence_id == "P1EN081" ] <- NA

# P2S90IDESUBUM		Ensete superbum
# P2S91IDEGLACUM	Ensete glaucum
# P2S94IDLECONG1	Ensete lecongkietii
# P2S95IDLECONG2	Ensete lecongkietii
# outgroup samples intially found on a long branch length, likely to bias phylogentic results
# type changed to NA 
md$type [ md$sequence_id == "P2S90IDESUBUM"  ] <- NA
md$type [ md$sequence_id == "P2S91IDEGLACUM" ] <- NA
md$type [ md$sequence_id == "P2S94IDLECONG1" ] <- NA
md$type [ md$sequence_id == "P2S95IDLECONG2" ] <- NA

# P1EN067 changed from Domesticated to Disease based on blast of bacterial genomes
md$type [ md$sequence_id == "P1EN067"  ] <- "Disease"

# P1EN068	Ferezye
# P1EN069	Nechwe
# changed from Domesticated to Disease
# wrongly classified accorinding to James
md$type [ md$sequence_id == "P1EN068"  ] <- "Disease"
md$type [ md$sequence_id == "P1EN069"  ] <- "Disease"

# remove Disease and NA samples
md <- md [ -which(md$type == "Disease" | is.na(md$type)) , ]

# leading zeros missing from sample_id
# sample_id also contains non-numeric id for some outgroup samples
# function to fix "sample_id"
# if non-numeric - NA
# if numeric - leading zeros added
fix_sample_id <- function(x) {
  x = suppressWarnings(as.numeric(x))
  if(!is.na(x)) {
    leading_zeros <- paste(rep("0", (6 - nchar(x))), collapse = "")
    return(paste0(leading_zeros, x))
  } else {
    return(NA)
  }
}

# see diff
# md$sample_id
# sapply(md$sample_id, fix_sample_id)

# change sample_id
md$sample_id <- unlist(sapply(md$sample_id, fix_sample_id))

# change Landrace_2 to landrace
colnames(md)[which(colnames(md) == "Landrace_2")] <- "landrace"

# correct some landrace names
md$landrace <- gsub("Eventricosum_Maurelii", "Maurelii (horticultural)", md$landrace)
md$landrace <- gsub("E.ventricosum seedling", "E.ventricosum", md$landrace)
md$landrace <- gsub("Elivingstonianum", "E.livingstonianum", md$landrace)
md$landrace <- gsub("ITC_sample", "ITC1387", md$landrace)
md$landrace <- gsub("white", "(white)", md$landrace)
md$landrace <- gsub("red", "(red)", md$landrace)
md$landrace <- gsub("green", "(green)", md$landrace)
md$landrace <- gsub("Wanadiye", "Wanade", md$landrace)
md$landrace [ startsWith(md$landrace, "Unknown") & md$type == "Semi-domesticated" ] <- "Unnamed (semi-domesticated)"
md$landrace [ startsWith(md$landrace, "Feral")   & md$type == "Semi-domesticated" ] <- "Unnamed (semi-domesticated)"
md$landrace [ startsWith(md$landrace, "Unknown") & md$type == "Domesticated" ]      <- "Unnamed"

# distinguish North and South Astara
# md [ md$landrace == "Astara" , ]
md$landrace [ md$sequence_id == "P1EN009" ] <- "Astara North"
md$landrace [ md$sequence_id == "P3S79ID602" ] <- "Astara South"
# md [ startsWith(md$landrace, "Astara") , ]

# change X_latitude to latitude
colnames(md)[which(colnames(md) == "X_latitude")] <- "latitude"

# change X_longitude to longitude
colnames(md)[which(colnames(md) == "X_longitude")] <- "longitude"

# drop Unique_code
md <- md[,-which(colnames(md) == "Unique_code")]

# change Sample_date to sample_date
colnames(md)[which(colnames(md) == "Sample_date")] <- "sample_date"

# get ethiopia regions
map_extract <- read.table("additional_data/map_extract.txt", header=T, sep = "\t")
map_extract <- select(map_extract, sequence_id, regions, subregions)

# remove commas 
map_extract$regions    <- gsub(",","", map_extract$regions)
map_extract$subregions <- gsub(",","", map_extract$subregions)

# translate N, S, E, W
map_extract$subregions <- gsub("Misraq","East",  map_extract$subregions)
map_extract$subregions <- gsub("Mirab" ,"West",  map_extract$subregions)
map_extract$subregions <- gsub("Debub" ,"South", map_extract$subregions)

# join map extract to md
md <- left_join(md, map_extract, by = "sequence_id")

# get index of wild
index_wild <- which(md$type == "Wild")
# md[index_wild,]

# add subregion to wild landrace name
md$landrace[index_wild] <- paste0("Wild (", md$subregions[index_wild], ")")
# md[index_wild,]

### add popmap data
popmap <- paste0("pop", 1:nrow(md))
md <- cbind(popmap, md)


### check sample number per type

# table(md$type)
# dim(md)

# head(md)

### write file
write.csv(md, "tGBS_metadata_phylogenetic_analysis.csv", row.names = F, quote=F)

