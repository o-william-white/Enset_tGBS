# run as
# Rscript plot_plink_pca.R <plink_path>

# use the function commandArgs to pass commands to the script
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

# read in pca data
eval <- read.table(file=paste0(args[1], "/plink.eigenval"), header=FALSE)
evec <- read.table(file=paste0(args[1], "/plink.eigenvec"), header=FALSE)

# set column names
colnames(evec) <- c("pop", "sample", paste0("e", 1:20))
colnames(eval) <- "eigenvalues"

# read in meta data with sample info
sample.metadata <- read.csv("/data/scratch/mpx469/tGBS_enset_project/tGBS_metadata_phylogenetic_analysis.csv")

# remove ourgroups
sample.metadata <- sample.metadata [ sample.metadata$type != "Outgroup" ,]

# change levels of type
sample.metadata$type <- factor(sample.metadata$type, levels=c("Domesticated", "Semi-domesticated", "Wild"))

# join dataframes to include all info
dat <- left_join(evec, sample.metadata, by=c("sample" = "sequence_id"))

# write data
write.table(dat, paste0(args[1], "/pca_dat.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t" )

# eigenvalues as percentages
# https://speciationgenomics.github.io/pca/
eval$pct <- eval$eigenvalues/sum(eval$eigenvalues)*100

# create palette for sample type
pal_sample <- brewer.pal(3, "Dark2")

# base plot
pdf(paste0(args[1], "/plot_pca.pdf"), height=10, width=10)
par(mfrow=c(2,2))

plot(dat$e1, dat$e2, col=dat$type,
   xlab=paste0("PC1 (", round(eval$pct[1],digits=2), " %)"),
   ylab=paste0("PC2 (", round(eval$pct[2],digits=2), " %)"))

legend("topright", inset=0.025, levels(dat$type), col=1:3, pch=21)


plot(dat$e1, dat$e3, col=dat$type,
   xlab=paste0("PC1 (", round(eval$pct[1],digits=2), " %)"),
   ylab=paste0("PC3 (", round(eval$pct[3],digits=2), " %)"))
   
plot(dat$e2, dat$e3, col=dat$type,
   xlab=paste0("PC2 (", round(eval$pct[2],digits=2), " %)"),
   ylab=paste0("PC3 (", round(eval$pct[3],digits=2), " %)"))

barplot(eval$eigenvalues, ylab="Eigenvalues", xlab="Principal components" )

dev.off()


png(paste0(args[1], "/plot_pca.png"), height=10, width=10, units="in", res=300)
par(mfrow=c(2,2))

plot(dat$e1, dat$e2, col=dat$type,
   xlab=paste0("PC1 (", round(eval$pct[1],digits=2), " %)"),
   ylab=paste0("PC2 (", round(eval$pct[2],digits=2), " %)"))

legend("topright", inset=0.025, levels(dat$type), col=1:3, pch=21)


plot(dat$e1, dat$e3, col=dat$type,
   xlab=paste0("PC1 (", round(eval$pct[1],digits=2), " %)"),
   ylab=paste0("PC3 (", round(eval$pct[3],digits=2), " %)"))
   
plot(dat$e2, dat$e3, col=dat$type,
   xlab=paste0("PC2 (", round(eval$pct[2],digits=2), " %)"),
   ylab=paste0("PC3 (", round(eval$pct[3],digits=2), " %)"))

barplot(eval$eigenvalues, ylab="Eigenvalues", xlab="Principal components" )

dev.off()



# ggplot

p1 <- ggplot(data=dat, aes(x=e1, y=e2, col=type)) +
   geom_point(alpha=0.5, size=2.5) +
   scale_color_manual(values=pal_sample) +
   labs(x=paste0("PC1 (", round(eval$pct[1],digits=2), " %)"), y=paste0("PC2 (", round(eval$pct[2],digits=2), " %)"), colour=NULL) +
   theme_classic() +
   theme(legend.position=c(0.9,0.9),
         legend.background = element_rect(colour="black"),
                 axis.text = element_text(colour="black"),
                 axis.line = element_line(colour="black"),
                 axis.ticks = element_line(colour="black"))
				 
p2 <- ggplot(data=dat, aes(x=e1, y=e3, col=type)) +
   geom_point(alpha=0.5, size=2.5) +
   scale_color_manual(values=pal_sample) +
   labs(x=paste0("PC1 (", round(eval$pct[1],digits=2), " %)"), y=paste0("PC3 (", round(eval$pct[3],digits=2), " %)"), colour=NULL) +
   theme_classic() +
   theme(legend.position="none",
         legend.background = element_rect(colour="black"),
                 axis.text = element_text(colour="black"),
                 axis.line = element_line(colour="black"),
                 axis.ticks = element_line(colour="black"))
				 
p3 <- ggplot(data=dat, aes(x=e2, y=e3, col=type)) +
   geom_point(alpha=0.5, size=2.5) +
   scale_color_manual(values=pal_sample) +
   labs(x=paste0("PC2 (", round(eval$pct[2],digits=2), " %)"), y=paste0("PC3 (", round(eval$pct[3],digits=2), " %)"), colour=NULL) +
   theme_classic() +
   theme(legend.position="none",
         legend.background = element_rect(colour="black"),
                 axis.text = element_text(colour="black"),
                 axis.line = element_line(colour="black"),
                 axis.ticks = element_line(colour="black"))

p4 <- ggplot(data=eval, aes(x=1:20,y=eigenvalues)) + geom_bar(stat="identity") +
   labs(x="Principal components", y="Eigenvalues") +
   theme_classic() +
   theme(legend.background = element_rect(colour="black"),
                 axis.text = element_text(colour="black"),
                 axis.line = element_line(colour="black"),
                 axis.ticks = element_line(colour="black"))		 


pdf(paste0(args[1], "/plot_pca_ggplot.pdf"), height=10, width=10)
plot_grid(p1,p2,p3,p4,nrow=2,ncol=2,align=TRUE)
dev.off()

png(paste0(args[1], "/plot_pca_ggplot.png"), height=10, width=10, units="in", res=300)
plot_grid(p1,p2,p3,p4,nrow=2,ncol=2,align=TRUE)
dev.off()


