
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(ape)
library(ggtree)
library(aplot)
library(RColorBrewer)

# read dsuite output
dsuite_BBAA <- read.table("dsuite/by_sample_BBAA.txt", header = TRUE)
dsuite_Dmin <- read.table("dsuite/by_sample_Dmin.txt", header = TRUE)
dsuite_tree <- read.table("dsuite/by_sample__tree.txt", header = TRUE)

# function to format and summarise the results
my_func <- function(dat) {
   
   # test
   dat = dsuite_BBAA
   
   # separate P2, P2 and P3 columns
   dat <- separate(dat, col = P1, into = c("P1_sequence_id", "P1_type"), sep = "_") %>% 
          separate(.,   col = P2, into = c("P2_sequence_id", "P2_type"), sep = "_") %>% 
          separate(.,   col = P3, into = c("P3_sequence_id", "P3_type"), sep = "_") 
   
   # check 
   # head(dat)
   # nrow(dat)

   # Benjamini-Hochberg (BH) correction for multiple tests
   dat$p.value.adj <- p.adjust(dat$p.value, method = "BH")
   
   # select P2 and P3 types for each test and sort
   tmp_type <- as.data.frame(t(apply(X = dat[,c("P2_type", "P3_type")], MARGIN = 1, FUN = sort)))
   colnames(tmp_type) <- c("Px", "Py")
   
   # get topology 
   tmp_type$topology_type <- paste(tmp_type$Px, tmp_type$Py)
   
   # add to dsuite table
   dat$topology_type <- tmp_type$topology_type
   
   # rm tmp
   rm(tmp_type)
   
   # check
   # head(dat)
   
   # select P2 and P3 sequence id for each test and sort
   tmp_sequence_id <- as.data.frame(t(apply(X = dat[,c("P2_sequence_id", "P3_sequence_id")], MARGIN = 1, FUN = sort)))
   colnames(tmp_sequence_id) <- c("Px", "Py")
   
   # get topology 
   tmp_sequence_id$topology_sequence_id <- paste(tmp_sequence_id$Px, tmp_sequence_id$Py)
   
   # add to dsuite table
   dat$topology_sequence_id <- tmp_sequence_id$topology_sequence_id
   
   # rm tmp
   rm(tmp_sequence_id)
   
   # as factor
   dat$topology_type <- as.factor(dat$topology_type)
   dat$topology_sequence_id <- as.factor(dat$topology_sequence_id)

   output <- vector(mode = "list", 4)
   names(output) <- c("total", "type", "sequence_id", "raw")
   
   output[[1]] <- dat %>% 
      summarise(number_of_tests             = n(),
                number_of_significant_tests = sum(p.value.adj < 0.01),
                pct_significant_tests       = round((sum(p.value.adj < 0.01) / n())*100, 2) )
   
   output[[2]] <- group_by(dat, topology_type) %>%
      summarise(number_of_tests             = n(),
                number_of_significant_tests = sum(p.value.adj < 0.01),
                pct_significant_tests       = round((sum(p.value.adj < 0.01) / n())*100, 2) )
   
   output[[3]] <- group_by(dat, topology_sequence_id, topology_type) %>%
      summarise(number_of_tests             = n(),
                number_of_significant_tests = sum(p.value.adj < 0.01),
                pct_significant_tests       = round((sum(p.value.adj < 0.01) / n())*100, 2),
                maxD = max(Dstatistic), 
                maxF = max(f4.ratio),
                minP = min(p.value.adj))
   
   output[[4]] <- dat
   
   output
   
}

# run for each dataset
results_BBAA <- my_func(dsuite_BBAA)
#results_Dmin <- my_func(dsuite_Dmin)
#results_tree <- my_func(dsuite_tree)

results_BBAA$total
results_BBAA$type
#results_BBAA$sequence_id
#results_BBAA$raw

# results_Dmin$total
# results_Dmin$type
# results_Dmin$sequence_id
# results_Dmin$raw
# 
# results_tree$total
# results_tree$type
# results_tree$sequence_id
# results_tree$raw

#nrow(dsuite_BBAA)
#nrow(results_BBAA$raw)

#head(dsuite_BBAA)
#head(results_BBAA$raw)

# write results
write.table(results_BBAA$type,        "dsuite/by_sample_dsuite_BBAA_type_summary.txt",        row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
write.table(results_BBAA$sequence_id, "dsuite/by_sample_dsuite_BBAA_sequence_id_summary.txt", row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
write.table(results_BBAA$raw,         "dsuite/by_sample_dsuite_BBAA_raw_summay.txt",          row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)

# get heatmap data
heatmap_dat <- separate(results_BBAA$type, topology_type, sep = " ", into = c("px", "py")) %>%
   select(px, py, pct_significant_tests)

# write heatmaps data
write.table(heatmap_dat, "dsuite/dsuitue_heatmap_data.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

# write png
png("dsuite/dsuite_summary_heatmap.png", height = 4.5, width = 4, units = "in", res = 600)
ggplot(heatmap_dat, aes(px, py, fill=pct_significant_tests)) + 
   geom_tile() + 
   scale_fill_gradient(low = "white", high = "red", na.value = "gray") + 
   geom_text(aes(label=paste(pct_significant_tests,"%"))) + 
   xlab("") +
   ylab("") +
   scale_y_discrete(expand = c(0,0)) +
   scale_x_discrete(expand = c(0,0)) +
   theme_bw() +
   theme(axis.line = element_line(colour = "black"),
         panel.background = element_rect(fill="lightgray", colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(), 
         axis.text.x = element_text(colour = "black", size=16, angle = 40, hjust = 1),
         axis.text.y = element_text(colour = "black", size=16),
         legend.position = "top") + 
   guides(fill = guide_colourbar(title = "Significant tests (%)",
                                 direction = "horizontal", 
                                 title.position = "top"))
dev.off()   

