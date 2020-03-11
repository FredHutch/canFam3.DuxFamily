library(DESeq2)
library(ensOrtholog.db)
data(mouse_human_ensOrtholog)

pkgDir <- "~/tapscott/RNA-Seq/canFam3.DuxFamily"
load(file.path(pkgDir, "data", "C2C12.ens.ddsl.rda")) ## HinM and MinM vs Luc
load(file.path(pkgDir, "data", "C2C12.ens.ddsl2.rda")) ## HinM and MinM vs NODOX
load(file.path(pkgDir, "data", "HinH.ens.dds.rda"))
     
#'
#' get an one-to-one mapping between mouse and human
#' 
load("~/tapscott/R_package/ensOrtholog.db/data/mouse_human_ensOrtholog.rda")
ortholog <- mouse_human_ensOrtholog

#' clean up and return one-to-one mapping
clean <- orthologMapping(ortholog=ortholog, domain.column="ensembl_gene_id",
                        image.column="hsapiens_homolog_ensembl_gene",
                        keep.id=NULL,
                        unique.image=TRUE)
i <- which(colnames(clean) == "domain:image")
colnames(clean)[i] <- "mouse.to.human.mapping"

#' mouse  DESeq results
ddsl <- c(C2C12.ens.ddsl, C2C12.ens.ddsl2)

C2C12.res <- lapply(names(ddsl), function(x) {
    res <- results(ddsl[[x]], lfcThreshold=1, alpha=0.1)
    colnames(res) <- paste0(x, ".", colnames(res))
    res
})
mouse <- do.call(cbind, C2C12.res)

#' human DESeq results
HinH.res <- results(HinH.ens.dds, lfcThreshold=1, alpha=0.1)
colnames(HinH.res) <- paste0("HinH.", colnames(HinH.res))

#' combine dataset
geneset <- intersect(rownames(mouse), rownames(clean))
sub.clean <- clean[geneset, ]
df <- append(as(sub.clean, "DataFrame"), mouse[geneset, ])

keep <- df$hsapiens_homolog_ensembl_gene %in% rownames(HinH.res)
df <- df[keep, ]
df <- append(df, HinH.res[df$hsapiens_homolog_ensembl_gene, ])

f <- file.path(pkgDir, "inst", "stats",
                             "mouse_human_homology_MasterSheet.csv")
write.csv(df, file=f)

#' how many multiple-to-multiple mapping?
i = which(!df$mouse.to.human.mapping=="1:1")
length(i) #449

#' example for Jenn
pkgDir <- "~/tapscott/RNA-Seq/canFam3.DuxFamily"
f <- file.path(pkgDir, "inst", "stats",
               "mouse_human_homology_MasterSheet.csv")
df <- read.csv(f, row.name=1)
df <- df[df$mouse.to.human.mapping=="1:1", ]
HinH <- df$HinH.log2FoldChange
HinM <- df$HinMvsLuc.log2FoldChange
MinM <- df$MinMvsLuc.log2FoldChange
getCor <- function(x, y, ...) {
    i <- is.na(x) | is.na(y)
    cor(x[!i], y[!i], ...)
}
getCor(HinH, HinM)
getCor(HinH, MinM)
getCor(MinM, HinM)
    
a=df$MinMvsLuc.log2FoldChange
b=df$MinMvsNODOX.log2FoldChange
getCor(a, b)

## The result contains mouse-human homology genes that have some expression
## in either mouse or human. This mean you might expect to see lfc=NA
## (no expression) in mosue or human, but not in both.

## The p-value is based upon the null hypthesis lfc = 1.

## Note that the original homology mapping is not one-to-one, many of them
## are one-to-multiple or multiple-to-multiple. Here, I simplied the mapping
## to 1(mouse)-to-1(human).If a mouse gene
## maps to multiple human genes, a human gene is selected by whichever comes
## first on the database and vice versa.


