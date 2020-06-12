
# run as
# Rscript bp_coverage_per_chr.R <blastn_output>

# use the function commandArgs to pass commands to the script
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(bedr)



# read blastn output
blast <- read.table(paste0("blastn_output/", args[1], "_blastn_out"), header=FALSE, sep="\t")

# set column names for blast output
colnames(blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue","bitscore", "staxid", "ssciname")



# read in all possible xcm contigs
chr <- read.table("contigs_xcm", col.names="chr", colClasses="character")

# create empty bp vector
bp <- vector(mode="numeric", length=length(chr))

# create df to populate later
df <- data.frame(chr, bp)




# if else loop to populate dataframe
if (filter(blast, ssciname=="Xanthomonas campestris pv. musacearum NCPPB 2005") %>% nrow == 0) {

        # if not blast hits to xcm do nothing, coverage values already set at 0

} else {

    # if there are reads identified as xcm, format a bed file
        
        bed <- filter(blast, ssciname=="Xanthomonas campestris pv. musacearum NCPPB 2005") %>% # filter for xcm hits
                  rowwise %>%                          # group rowise
                  mutate(start = min(sstart, send),    # find min from sstart and send (blast can be in reverse complement)
                    end = max(sstart, send)) %>%       # find max from sstart and send
                  ungroup() %>%                        # ungroup to work colwise again
                  select(sseqid, start, end) %>%       # select columns for sseqid, start, send
                  rename("chr" = "sseqid") %>%         # rename "sseqid" as "chr"
                  mutate(chr = as.character(chr)) %>%  # change chr to character class
                  distinct %>%                         # reomve duplicate rows
                  arrange(chr, start)                  # arrange
        
        # write bed file if present
        write.table(bed, paste0("summary_table_output/xcm_bed/", args[1], ".bed"), col.names=FALSE, quote=FALSE, sep="\t", row.names=FALSE)
        
    # for each posibe chromosome in xcm run bedr to merge blast hits and sum coverage
        for (i in df$chr) {

           # if chr not present in bed file do nothing, coverage values already set at 0

           if (filter(bed, chr == i) %>% nrow == 0) {

           } else {

           # if chr present, merge overlapping hits and count bp coverage
           
           bed.i <- filter(bed, chr == i)

           bp.i <- bedr.merge.region(data.frame(bed.i), check.chr = FALSE, verbose=FALSE) %>%    # bedr to merge overlapping blast hits
                      mutate(length = end - start) %>%                                                       # calculate length of each region
                      summarise(sum = sum(length)) %>%                                                       # sum lengths
                      as.numeric
           
           # populate dataframe
           df[ df$chr == i , "bp" ] <- bp.i

           }
    }
}

write.table(df, paste0("summary_table_output/xcm_cov/", args[1], "_cov"), col.names=TRUE, quote=FALSE, sep="\t", row.names=FALSE)
