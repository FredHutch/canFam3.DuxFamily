# This script performs repeat element analysis

pkg_dir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/canFam3.DuxFamily"
rmsk_dir <- "/fh/fast/tapscott_s/CompBio/R_packages/rmskStats"
fig_dir <- file.path(pkg_dir, "figures")
data_dir <- file.path(pkg_dir, "data")

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)

load(file.path(pkg_dir, "data", "rmsk.dds.rda"))
load(file.path(pkg_dir, "data", "rmsk.res.rda"))

#
# DEseq2 resutls + MAplot + scatter plot
#
rmsk_HinC <- rmsk.dds[["HinC"]]
rmsk_CinC <- rmsk.dds[["CinC"]]
rmsk_HinC_df <- as(results(rmsk_HinC, alpha=0.05, lfcThreshold = 0.5), "data.frame")
rmsk_CinC_df <- as(results(rmsk_CinC, alpha=0.05, lfcThreshold = 0.5), "data.frame")

# scatter plot


# What threshold to choose. If I choose different threshold and Null hypothesis, 
# will I got the some answer?

# HinC
as(results(rmsk_HinC, alpha=0.05), "data.frame") %>%
  summarise(up = sum(padj < 0.05 & log2FoldChange > 0.95, na.rm=TRUE), 
            down = sum(padj < 0.05 & log2FoldChange < 0.95, na.rm=TRUE))

as(results(rmsk_CinC, alpha=0.05), "data.frame") %>%
  summarise(up = sum(padj < 0.05 & log2FoldChange > 0.95, na.rm=TRUE), 
            down = sum(padj < 0.05 & log2FoldChange < 0.95, na.rm=TRUE))


res <- results(rmsk_HinC, alpha=0.05, lfcThreshold = 0.5)
summary(res)
pdf(file.path(fig_dir, "rmsk_HinC_MAplot.pdf"), height=4, width=4.5)
plotMA(res, main="HinC: 726 annotated repeat elements")
text(x=10000, y=4, labels="up: 72", adj=c(0, 0))
text(x=10000, y=-4, labels="down: 46", adj=c(0, 1))
dev.off()

# CinC
res <- results(rmsk_CinC, alpha=0.05), lfcThreshold = 0.5)
summary(res)

pdf(file.path(fig_dir, "rmsk_CinC_MAplot.pdf"), height=4, width=4.5)
plotMA(res, main="CinC: 717 annotated repeat elements")
text(x=10000, y=2.5, labels="up: 58", adj=c(0, 0))
text(x=10000, y=-2.5, labels="down: 26", adj=c(0, 1))
dev.off()