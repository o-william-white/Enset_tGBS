
### load libraries
library(vcfR)
library(adegenet)
library(reshape2)
library(ggplot2)
library(dplyr)
library(cowplot)
library(poppr)

args = commandArgs(trailingOnly=TRUE)


### read vcf

# read vcf file
vcf <- read.vcfR(paste0("populations/populations_select_80_all_snps_whitelist_blacklist_", args[1], "/populations.snps.vcf"))

# get alt sites
site_alt <- paste(getID(vcf), getALT(vcf), sep=".")



### read sample metadata

# read data
sample_metadata <- read.csv("tGBS_metadata_population_genetics.csv", colClasses = "character")

# check numbers per type
# table(sample_metadata$type)

# check identical
print("Checking vcf and sample metadata names are identical")
identical(colnames(vcf@gt)[-1], sample_metadata$sequence_id)


### get allele annotations 

# read annotations data in wide format
annotations <- read.table(paste0("snpeff/", args[1], "/snpsift_wide_ann_fmt.txt"), header = TRUE, sep="\t")

# nrow(annotations)

# note the annotations are for the alt allele and their putative effect
# head(annotations)
# head(site_alt)

# get all unique annotations to order levels
annotations_levels <- sort(unique(melt(annotations, id.vars="Site")$value))



### create genind object with population info

# convert vcf to genind
gi <- vcfR2genind(vcf, return.alleles=TRUE)

# check identical
print("Checking sample metadata and genind names are identical")
identical(sample_metadata$sequence_id, indNames(gi))

# add population info
gi@pop <- as.factor(sample_metadata$type)

# check numbers per type
# table(gi$pop)



### get samples for each population

# get domesticated and wild samples
samples_domesticated <- indNames(popsub(gi, "Domesticated"))
samples_wild         <- indNames(popsub(gi, "Wild"))



### create empty lists to populate

# counts allele categories
list_count <- vector("list", 100)

# retained
list_retained <- vector("list", 100)

# novel
list_novel <- vector("list", 100)

# counts annotations
list_annotations_shared               <- vector("list", 100)
list_annotations_private_domesticated <- vector("list", 100)
list_annotations_private_wild         <- vector("list", 100)

# proportion annotations
list_annotations_shared_proportion               <- vector("list", 100)
list_annotations_private_domesticated_proportion <- vector("list", 100)
list_annotations_private_wild_proportion         <- vector("list", 100)


for(i in 1:100) {
  
  # testing
  #i=1
  
  # sample 13 random domesticated samples - equal in length to the wild samples
  samples_domesticated_i <- sample(samples_domesticated, length(samples_wild), replace = FALSE)
  
  # subset gi for domesticated sample selection
  gi_domestic_i <- gi[ match(samples_domesticated_i, rownames(gi$tab)) ,  ]
  
  # create df of allele presence in domesticated and wild
  df <- data.frame(Domesticated = apply(gi_domestic_i$tab, 2, function(x) any(x!=0, na.rm = TRUE)), 
                   Wild         = apply(popsub(gi, "Wild", drop = FALSE)$tab, 2, function(x) any(x!=0, na.rm = TRUE)))
  
  # move rownames to column
  df <- cbind(Site = row.names(df), df)
  row.names(df) <- NULL

  # subset df for alt sites
  df <- df [ match(site_alt, df$Site) , ]
  
  # check
  # head(df)

  ### categorise alleles
  
  # categorised as:
  #    - shared between domesticated and wild
  #    - private to domesticated
  #    - private to wild
  df$category <- ifelse(df$Domesticated == TRUE  & df$Wild == TRUE,  "Shared", 
                 ifelse(df$Domesticated == TRUE  & df$Wild == FALSE, "Private domesticated", 
                 ifelse(df$Domesticated == FALSE & df$Wild == TRUE,  "Private wild", 
                 ifelse(df$Domesticated == FALSE & df$Wild == FALSE,  "Missing", NA))))

  # category as factor
  df$category <- factor(df$category, levels = c("Shared", "Private domesticated", "Private wild", "Missing"))

  # check 
  # head(df)



  ### counts per category

  # counts
  summary_allele_categories <- table(df$category)

  list_count[[i]] <- summary_allele_categories

  
  
  ### proportion of wild diversity retained 

  # calculated as (shared / total found in wild)
  retained <- unname(summary_allele_categories["Shared"] / (summary_allele_categories["Shared"] + summary_allele_categories["Private wild"]))
  
  names(retained) <- "Retained"
  
  list_retained[[i]] <- retained


  
  ### novel domeseticated diversity

  # calculated as (domesticted / total)
  novel <- summary_allele_categories["Private domesticated"] /  sum(summary_allele_categories["Shared"] + summary_allele_categories["Private domesticated"] + summary_allele_categories["Private wild"])

  names(novel) <- "Novel"

  list_novel[[i]] <- novel



  ### subset df for alt alleles (only alt alleles are annotated with effects) and join annotations
  
  # subset df for alt sites
  # df <- df [ match(site_alt, df$Site) , ]
  
  # join annotations
  df <- left_join(df, annotations, by = "Site")
  
  # check 
  # head(df)



  ### get counts and proportion of annotation types for alleles in each category
  
  summarise_ann <- function(x) {
    tmp <- filter(df, category == x) %>%                          # filter for category
      select(., Site, starts_with("EFF")) %>%                         # select sites and effects columns
      melt(id.vars="Site") %>%                                        # melt to long
      mutate(value = factor(value, levels = annotations_levels))      # set levels for annotations
    table(tmp$value)
  }

  # counts
  list_annotations_shared[[i]]               <- summarise_ann("Shared")
  list_annotations_private_domesticated[[i]] <- summarise_ann("Private domesticated")
  list_annotations_private_wild[[i]]         <- summarise_ann("Private wild")

  # proportion
  list_annotations_shared_proportion[[i]]               <- summarise_ann("Shared")               / summary_allele_categories["Shared"]
  list_annotations_private_domesticated_proportion[[i]] <- summarise_ann("Private domesticated") / summary_allele_categories["Private domesticated"]
  list_annotations_private_wild_proportion[[i]]         <- summarise_ann("Private wild")         / summary_allele_categories["Private wild"]

  print(paste("rep", i, "complete"))
  
}



# function to summarise lists
summarise_list <- function(x) {
  # create matrix for each
  m <- do.call("rbind", x)
  # produce df
  data.frame(name = colnames(m),
             mean = unname(colMeans(m)),
             sd   = unname(apply(m, 2, sd)))
}


# summarise lists
rf_count <- summarise_list(list_count)

rf_retained <- summarise_list(list_retained)

rf_novel <- summarise_list(list_novel)

rf_annotations_shared               <- summarise_list(list_annotations_shared)
rf_annotations_private_domesticated <- summarise_list(list_annotations_private_domesticated)
rf_annotations_private_wild         <- summarise_list(list_annotations_private_wild)

rf_annotations_shared_proportion               <- summarise_list(list_annotations_shared_proportion) 
rf_annotations_private_domesticated_proportion <- summarise_list(list_annotations_private_domesticated_proportion)
rf_annotations_private_wild_proportion         <- summarise_list(list_annotations_private_wild_proportion)

pal <- c("#0072B2", "#D55E00", "#E69F00", "black")
names(pal) <- c("Domesticated", "Feral", "Wild", "Outgroup")
#pal

pal_private <- c("gray", "#0072B2",  "#E69F00")
names(pal_private) <- c("Shared", "Private domesticated", "private wild")
#pal_private


### plot mean counts based on rarefaction

# drop missing
rf_count <- rf_count[-4,]

# allele category as factor
rf_count$name <- factor(rf_count$name, levels = c("Shared", "Private domesticated", "Private wild" ))

# plot
plot_rf_count <- ggplot(rf_count, aes(x=name, y = mean)) +
  geom_bar(stat="identity", fill=pal_private) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1*max(rf_count$mean))) + 
  xlab("") + 
  ylab("Alleles (mean N)") + 
  theme_classic() + 
  theme(axis.text = element_text(colour = "black"))  

# check 
#plot_rf_count



### plot annotation effects based on rarefaction

# function to format annotation data from rf
fmt_annotation_rf <- function(shared, private_domesticated, private_wild) {
  
  # rbind to long df
  tmp <- rbind(shared,private_domesticated,private_wild)
  
  # add column for allele category
  tmp$category <- c(rep("Shared",               nrow(shared)),
                    rep("Private domesticated", nrow(private_domesticated)),
                    rep("Private wild",         nrow(private_wild)))
  
  # set levels for annotations
  tmp$name <- factor(tmp$name, levels = annotations_levels)
  
  # set levels for allele category
  tmp$category <- factor(tmp$category, levels = c("Shared", "Private domesticated", "Private wild" ))
  
  # edit annotations names for plotting
  #   - change "_" to " "
  #   - change "&" to " & "
  tmp$name <- gsub("_", " ",   tmp$name)
  tmp$name <- gsub("&", " & ", tmp$name)

  # return
  tmp  
}

# run function
rf_annotations_all <- fmt_annotation_rf(rf_annotations_shared,
                                        rf_annotations_private_domesticated,
                                        rf_annotations_private_wild)

rf_annotations_all_proportion <- fmt_annotation_rf(rf_annotations_shared_proportion,
                                                   rf_annotations_private_domesticated_proportion,
                                                   rf_annotations_private_wild_proportion)

# plot
plot_rf_annotations <- ggplot(data=rf_annotations_all, aes(x=name, y=mean, fill=category)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, position = position_dodge(0.9)) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_grey() +
  xlab("") + 
  ylab("Mean number of annotations") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 315, hjust = 0)) +
  guides(fill = guide_legend(title = "Allele category"))

# check
#plot_rf_annotations


# plot proportion 
plot_rf_annotations_proportion <- ggplot(data=rf_annotations_all_proportion, aes(x=name, y=mean, fill=category)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, position = position_dodge(0.9)) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_grey() +
  xlab("") + 
  ylab("Mean number of annotations") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 315, hjust = 0)) +
  guides(fill = guide_legend(title = "Allele category"))

#plot_rf_annotations_proportion


# repeat for proportion of missense and synonymous only  

# subset df
rf_annotations_missense_proportion   <- rf_annotations_all_proportion [ rf_annotations_all_proportion$name == "missense variant"   , ]
rf_annotations_synonymous_proportion <- rf_annotations_all_proportion [ rf_annotations_all_proportion$name == "synonymous variant" , ]

# plot
plot_rf_annotations_missense_proportion <- ggplot(data=rf_annotations_missense_proportion, aes(x=category, y=mean)) + 
  geom_bar(stat="identity", fill=pal_private) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.2*max(rf_annotations_missense_proportion$mean))) + 
  xlab("") + 
  ylab("Proportion of missense annotations") + 
  theme_classic() + 
  theme(axis.text = element_text(colour = "black"))

plot_rf_annotations_synonymous_proportion <- ggplot(data=rf_annotations_synonymous_proportion, aes(x=category, y=mean)) + 
  geom_bar(stat="identity", fill=pal_private) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.2*max(rf_annotations_synonymous_proportion$mean))) + 
  xlab("") + 
  ylab("Proportion of synonymous annotations") + 
  theme_classic() + 
  theme(axis.text = element_text(colour = "black"))

plot_overall <- plot_grid(
  plot_rf_count + scale_x_discrete(labels = c("Shared", "Private\ndomesticated", "Private\nwild")) ,
  plot_rf_annotations_missense_proportion + scale_x_discrete(labels = c("Shared", "Private\ndomesticated", "Private\nwild")),
  align = "hv", nrow = 1)


# create output dir
dir.create(paste0("population_genetics/", args[1]), recursive = T, showWarnings = F)

# plot
png(paste0("population_genetics/", args[1], "/allele_categories.png"), height = 3, width = 6, res = 600, units = "in")
plot_overall
dev.off()

# write data
write.table(rf_count,                           paste0("population_genetics/", args[1], "/allele_categories_counts.txt"),   col.names = F, row.names = F, quote = F, sep = "\t")
write.table(rf_annotations_missense_proportion, paste0("population_genetics/", args[1], "/allele_categories_missense.txt"), col.names = F, row.names = F, quote = F, sep = "\t")

write.table(rf_retained,                        paste0("population_genetics/", args[1], "/allele_categories_retained.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(rf_novel,                           paste0("population_genetics/", args[1], "/allele_categories_novel.txt"),    col.names = F, row.names = F, quote = F, sep = "\t")

