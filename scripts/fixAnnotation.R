#'
#' (1) Add annotation from biomart (including homology)
#'
library(TxDb.Cfamiliaris.UCSC.canFam3.ensGene)
data(keys) ## return keys
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/canFam3.DuxFamily"
csvFiles <- list.files(file.path(pkgDir, "inst", "stats"), pattern="dds.csv",
                       full.names=TRUE)
csvFiles <- csvFiles[grep("ens_", csvFiles)]
lapply(csvFiles, function(x) {
    df <- read.csv(x, row.names=1)
    anno <- keys[rownames(df), ]
    df$SYMBOL <- anno$external_gene_name
    df$ENTREZID <- anno$entrezgene
    df$DESCRIPTION <- anno$description
    df <- cbind(df, anno[, 7:10])
    df <- df[, -4]
    df <- df[, c(1:9, 17:20, 10:16)]
    filename <- file.path(dirname(x),
                          paste0(sub(".csv", "", basename(x)), "_ensembl87.csv"))
    write.csv(df, file=filename)
})


