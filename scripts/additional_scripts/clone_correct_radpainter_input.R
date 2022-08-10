
# read radpainter input
rp <- read.table("radpainter/populations.haps.radpainter", header = T, sep = "\t")

# read mlg input
mlg <- read.table("estimate_mlg/mlg_farthest_bitwise_monophyletic_single_rep.txt", header = T, sep = "\t")

# filter outgroups from mlg
mlg <- mlg [ mlg$type != "Outgroup" , ]

# filter rp for single representatives per mlg
rp <- rp [ , match(mlg$sequence_id, colnames(rp)) ]

# colnames of rp match order of mlg
identical(as.character(colnames(rp)), as.character(mlg$sequence_id))

# write
write.table(rp, "radpainter/populations.haps.clone.correct.radpainter", col.names = T, row.names = F, quote = F, sep = "\t")

