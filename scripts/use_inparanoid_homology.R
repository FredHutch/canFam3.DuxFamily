#' This script use Bioconductor's homology package built by inparanoid algorithm
#' We need human to canine

library(DBI)
library("hom.Hs.inp.db")
library(dplyr)
library(tidyr)
library(purrr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Cf.eg.db)

# "hom.Hs.inpCANFA: human to canine
pkg_dir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/canFam3.DuxFamily"
source(file.path(pkg_dir, "scripts", "inparanoid_homology.R"))
load(file.path(pkg_dir, "data", "C2C12.ens.ddsl2.rda"))
load(file.path(pkg_dir, "data", "CinC.ens.dds.rda"))
load(file.path(pkg_dir, "data", "HinC.ens.dds.rda"))
load(file.path(pkg_dir, "data", "HinH.ens.dds.rda"))
load(file.path(pkg_dir, "data", "human_cleavage.rda"))


#
# cleavage_homology and 2C-lie
#
cleavage_homology <- human_inparanoid_homology(as.character(human_cleavage$GeneID))
# how many of them are differentially expressed? run GSEA
#like2c_homology

#
# whole gene set homology
#
load(file.path(pkg_dir, "data", "HinH.ens.SE.rda"))
no_keep <- grepl("LRG", rownames(HinH.ens.SE))
human_universe <- rownames(HinH.ens.SE)[!no_keep]
universe_homology <- human_inparanoid_homology(human_universe)

#' 05_compare_homology: compare homologues expression and functions
# (a) CinC and HinH expression on human homology
# (b) MinM and HinH expression on human homology
# (c) CinC and HinC expression on human homology
# use (a) and (b) to compare correlation 
library(DESeq2)
library(tidyverse)
HinH_res <- as(results(HinH.ens.dds, lfcThreshold=1, alpha=0.05), "data.frame") %>%
  rownames_to_column(var="HUMAN_ENSEMBL") %>%
  dplyr::select(HUMAN_ENSEMBL, log2FoldChange, padj) %>%
  rename(HinH_logFC=log2FoldChange, HinH_padj=padj)
CinC_res <- results(CinC.ens.dds, lfcThreshold=1, alpha=0.05)
HinC_res <- results(HinC.ens.dds, lfcThreshold=1, alpha=0.05)
MinM_res <- results(C2C12.ens.ddsl2[[2]], lfcThreshold=1, alpha=0.05)

HinH_CinC_HinC <- universe_homology[["human2canine"]] %>%
  left_join()

# combine human2canine and human2mouse
HinH_MinM <- universe_homology[[2]] %>%


# (a) build a matrix of LogFC
# (b) exame pairwise correlation

#
# cleavage-stage gene signature only
#



