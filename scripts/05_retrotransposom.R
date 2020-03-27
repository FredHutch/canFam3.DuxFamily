# This script performs repeat element analysis

pkg_dir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/canFam3.DuxFamily"
rmsk_dir <- "/fh/fast/tapscott_s/CompBio/R_packages/rmskStats"
fig_dir <- file.path(pkg_dir, "figures")
data_dir <- file.path(pkg_dir, "data")

library(DESeq2)
library(tidyverse)
library(tidyr)


load(file.path(pkg_dir, "data", "rmsk.dds.rda"))
load(file.path(pkg_dir, "data", "rmsk.res.rda"))

# summary of repFamily and 
# 
# DESeq2 results: based on two thresholds
#
repeat_hypergeometric <- function(universe, selected, x, k) {
    # x: total DE repName in the set (repClass or repFamily)
    # k: total repName in the set (repClass or repFamily)
    m <- selected # white ball
    n <- universe - selected # black ball
    prob <- dhyper(x=x, m=m, n=n, k=k)
    mu <- as.numeric(k * (m/(m+n)))
    return(c(prob=prob, mu=mu))
}

res_df <- lapply(names(rmsk.dds), function(name) {
    dds <- rmsk.dds[[name]]
    rowdata <- as(rowData(dds), "data.frame") %>%
      rownames_to_column(var="repName") %>%
      dplyr::select(repName, repClass, repFamily)
    # first threshold: p-value=0.0r correspondng to H_0: |lfc| < 0.5  
    res_df <- as(results(dds, alpha=0.05, lfcThreshold=0.5), "data.frame") %>%
      rownames_to_column(var="repName") %>%
      left_join(rowdata, by="repName") 
    # padj: adj p-value =0.05 corresponding to H_0: |lfc| < 0.5
    # padj_null_zero: adj p-value=0.05 correspondng to H_0: lfc=0  
    universe <- nrow(res_df)
    selected_up <- res_df %>% 
      summarise(up = sum(padj < 0.05 & log2FoldChange > 0.5, na.rm=TRUE)) %>%
      pull(up)
    selected_down <- res_df %>% 
      summarise(down = sum(padj < 0.05 & log2FoldChange < 0.5, na.rm=TRUE)) %>%
      pull(down)

    # repFamily enrichment/depletion summary
    repFamily_summary <- res_df %>% group_by(repFamily) %>%
      summarise(total=n(), 
                up = sum(padj < 0.05 & log2FoldChange > 0, na.rm=TRUE), 
                down = sum(padj < 0.05 & log2FoldChange < 0, na.rm=TRUE),
                enrichment_prob = dhyper(x=up, k=total, m=selected_up, n=universe - selected_up),
                enrichment_mu = total * selected_up/universe,
                depletion_prob = dhyper(x=down, k=total, m=selected_down, n=universe-selected_down),
                depletion_mu = total * selected_down/universe)

    # repClass enrichement/depletion summary
    repClass_summary <- res_df %>% group_by(repClass) %>%
       summarise(total=n(), 
                up = sum(padj < 0.05 & log2FoldChange > 0, na.rm=TRUE), 
                down = sum(padj < 0.05 & log2FoldChange < 0, na.rm=TRUE),
                enrichment_prob = dhyper(x=up, k=total, m=selected_up, n=universe - selected_up),
                enrichment_mu = total * selected_up/universe,
                depletion_prob = dhyper(x=down, k=total, m=selected_down, n=universe-selected_down),
                depletion_mu = total * selected_down/universe)      
    # print out to a spread sheet                
    return(list(res_df=res_df, repFamily_summary=repFamily_summary, repClass_summary=repClass_summary))               
})
names(res_df) <- names(rmsk.dds)



#
# (a) scatter of log2FoldChange in HinC and CinC, which 
#


# (a) define differentially expressed element

# 1. p-value = 0.05 corresponding to H_0: |lfc| < 0.5
# 2. p-value = 0.05 corresponding to H_0: lfc = 0 and cutoff at lfc = 0.95

# What threshold to choose? If I choose different threshold and Null hypothesis, 
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