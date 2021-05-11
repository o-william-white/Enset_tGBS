
# run as
# Rscript identify_duplicate_loci.R vcf output

# use the function commandArgs to pass commands to the script
args = commandArgs(trailingOnly=TRUE)

# get args
vcf    <- args[1]
output <- args[2]

library(vcfR)

# read vcf
vcf <- read.vcfR(file = vcf)

# create output dir
dir.create(output)

# dataframe with vcf chr, pos and id
df <- data.frame(chr = getCHROM(vcf),
                 pos = getPOS(vcf), 
                 id = getID(vcf))

# get ids for sites with identical chr and pos
duplicated_id <- df$id [ duplicated(paste0(df$chr, df$pos, sep="_")) ]

# get loci from duplicated id
duplicated_loci <- sapply(strsplit(as.character(duplicated_id), ":"), getElement, 1)

# write duplicated loci
write.table(duplicated_loci, paste0(output,"/duplicate_loci_blacklist.txt"), col.names = F, row.names = F, quote = F)

