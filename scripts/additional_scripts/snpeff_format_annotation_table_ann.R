
library(dplyr)
library(tidyr)
library(ggplot2)


### long data

# read long data
snpsift_long <- read.table("snpsift_long_ann.txt", header = TRUE, sep = "\t")

# check 
# head(snpsift_long)

# fmt 
snpsift_long <- snpsift_long %>% 
  rename("ANN" = "ANN....EFFECT") %>%
  mutate(Site = paste(ID, ALT, sep=".")) %>%
  select(Site, ANN)

# check again
# head(snpsift_long)


### wide format

# read wide data
snpsift_wide <- read.table("snpsift_wide_ann.txt", header = TRUE, sep = "\t")

# check 
# head(snpsift_wide)

# max number of annotations
# max(table(snpsift_long$Site))

# fmt
snpsift_wide <- snpsift_wide %>% 
  rename("ANN" = "ANN....EFFECT") %>%
  mutate(Site = paste(ID, ALT, sep=".")) %>% 
  select(Site, ANN) %>% separate(ANN, into = paste("EFF", 1:max(table(snpsift_long$Site)), sep = "_"), sep=",") 

# check 
# head(snpsift_wide)

# same number of sites
print("checking the length of snps in long and wide dataset are the same")
identical(length(unique(snpsift_long$Site)), length(unique(snpsift_wide$Site)))


### write formated files

write.table(snpsift_long, "snpsift_long_ann_fmt.txt", col.names = TRUE, row.names = FALSE, sep="\t", quote = FALSE)  
write.table(snpsift_wide, "snpsift_wide_ann_fmt.txt", col.names = TRUE, row.names = FALSE, sep="\t", quote = FALSE)  


### plot annotation counts

# see count for each type
# table(snpsift_long$ANN)

# as data.frame
df_eff <- as.data.frame(table(snpsift_long$ANN))

# sort levels
df_eff$Var1 <- factor(df_eff$Var1, levels = df_eff$Var1[ order(df_eff$Freq, decreasing = TRUE) ])

# plot frequency per annotation
p <- ggplot(data = df_eff, aes(x=Var1, y=Freq)) + geom_bar(stat = "identity") + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_classic() + 
  xlab("Annotation") + 
  ylab("Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"))

# save plot
png("plot_snpeff_ann.png", height=8, width = 12, units = "in", res = 600)
p
dev.off()
