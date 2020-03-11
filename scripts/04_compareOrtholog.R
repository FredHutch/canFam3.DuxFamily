library(DESeq2)
pkgDir <- "~/tapscott/RNA-Seq/canFam3.DuxFamily"
dataDir <- file.path(pkgDir, "data")
statsDir <- file.path(pkgDir, "inst", "stats")
library(TxDb.Cfamiliaris.UCSC.canFam3.ensGene)
data(keys)
     
#' Create master homology sheet with canine, mouse and human. The DESeq2 results is based
#' uopn the null hypothesis in which lfc > 1 and alpha = 0.1 for all canine, mouse and
#' human. The genes of canine that have some expression ( > 6 across samples) are
#' subject to DESeq analysis. So as human and mouse.

#'
#' load canine SE, dds and load human and mouse dds
#'
load(file.path(dataDir, "C2C12.ens.ddsl.rda"))
load(file.path(dataDir, "HinH.ens.dds.rda"))
load(file.path(dataDir, "canine.ens.SE.rda"))
## k9 dds
se <- canine.ens.SE[rowSums(assay(canine.ens.SE)) > 6, 1:9]
k9.ens.ddsl <- lapply(c("HinC", "CinC"), function(x) {
    sub <- se[, se$nick_name%in% c(x, "LinC")]
    print(dim(sub))
    sub$nick_name <- factor(sub$nick_name, levels=c("LinC", x))
    print(levels(sub$nick_name))
    dds <- DESeqDataSet(sub, design = ~ nick_name)
    dds <- DESeq(dds, parallel=TRUE)
})
names(k9.ens.ddsl) <- c("HinC", "CinC")

#'
#' Master sheet
#' 
k9.res <- lapply(names(k9.ens.ddsl), function(x) {
    res <- results(k9.ens.ddsl[[x]], lfcThreshold=1, alpha=0.1)
    colnames(res) <- paste0(x, ".", colnames(res))
    res
})
k9.keys <- keys[rownames(k9.res[[1]]), ]
k9 <- append(do.call(cbind, k9.res), as(k9.keys, "DataFrame"))
write.csv(k9, file=file.path(statsDir, "canine_MasterSheet.csv"))

#' k9 and mouse
names(C2C12.ens.ddsl) <- c("HinM", "MinM")
C2C12.res <- lapply(names(C2C12.ens.ddsl), function(x) {
    res <- results(C2C12.ens.ddsl[[x]], lfcThreshold=1, alpha=0.1)
    colnames(res) <- paste0(x, ".", colnames(res))
    res
})
mouse <- do.call(cbind, C2C12.res)

#' k9 and mouse
idx <- !k9$mmusculus_homolog_ensembl_gene == ""
df <- k9[idx, c("HinC.log2FoldChange", "HinC.padj",
                "CinC.log2FoldChange", "CinC.padj",
                "external_gene_name",
                "mmusculus_homolog_ensembl_gene",
                "mmusculus_homolog_associated_gene_name")]
df <- append(df,  mouse[df$mmusculus_homolog_ensembl_gene, c(2, 6, 8, 11)])
write.csv(df, file=file.path(statsDir, "canine_mouse_homology_MasterSheet.csv"))

#' k9 and human
HinH.res <- results(HinH.ens.dds, lfcThreshold=1, alpha=0.1)
colnames(HinH.res) <- paste0("HinH.", colnames(HinH.res))

idx <- !k9$hsapiens_homolog_ensembl_gene == ""
df <- k9[idx, c("HinC.log2FoldChange", "HinC.padj",
                "CinC.log2FoldChange", "CinC.padj",
                "external_gene_name",
                "hsapiens_homolog_ensembl_gene",
                "hsapiens_homolog_associated_gene_name")]
df <- append(df, HinH.res[df$hsapiens_homolog_ensembl_gene, c(2, 6)])
write.csv(df, file=file.path(statsDir, "canine_human_homology_MasterSheet.csv"))         

#' k9 and human and mouse
idx <- !k9$hsapiens_homolog_ensembl_gene == "" & !k9$mmusculus_homolog_ensembl_gene == ""
df <- k9[idx, c("HinC.log2FoldChange", "HinC.padj",
                "CinC.log2FoldChange", "CinC.padj",
                "external_gene_name",
                "hsapiens_homolog_ensembl_gene",
                "hsapiens_homolog_associated_gene_name", 
                "mmusculus_homolog_ensembl_gene",
                "mmusculus_homolog_associated_gene_name")]
df <- append(df, HinH.res[df$hsapiens_homolog_ensembl_gene, c(2, 6)])
df <- append(df, mouse[df$mmusculus_homolog_ensembl_gene, c(2, 6, 8, 11)])
write.csv(df, file=file.path(statsDir, "canine_human_mouse_MasterSheet.csv"))
             

## human and mouse (later)
