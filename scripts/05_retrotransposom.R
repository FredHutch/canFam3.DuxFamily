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

# scatter plot of repeat family counts
# scatter plot of repeat class 