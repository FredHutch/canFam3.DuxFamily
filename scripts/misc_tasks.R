# This script does misc. tasks
pkg_dir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/canFam3.DuxFamily"
libary(xlsx)

#'
#' Task 1: extract  cleavage-stage genes from Cairn's paper
#'
cairn_gene <- read.xlsx(file.path(pkg_dir, "inst", "extdata"))

#'
#' Task 2: extract mouse 2C-like genes from Akiyami paper
#'