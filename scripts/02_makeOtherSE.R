library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)
library(DESeq2)
bpparam=MulticoreParam(worker=6)

#'
#' C2C12 mm10: 
#'
pkgDir <- "~/tapscott/RNA-Seq/canFam3.DuxFamily"
dataDir <- file.path(pkgDir, "data")

ngsDir <- "/shared/ngs/illumina/jwhiddon/150918_SN367_0553_AH7YLFBCXX"
bamDir <- file.path(ngsDir, "tophat", "bam")
bamFiles <- list.files(bamDir, pattern="\\.bam$", include.dirs=TRUE,
                       all.files=FALSE, full.names=TRUE)
#' count hits
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
features <- exonsBy(TxDb.Mmusculus.UCSC.mm10.ensGene, by="gene")
features <- keepStandardChromosomes(features, pruning.mode="fine")
C2C12.ens.SE <- summarizeOverlaps(features=features,
                                  reads=bamFiles,
                                  mode="IntersectionStrict",
                                  inter.feature=TRUE,
                                  singleEnd=TRUE,
                                  ignore.strand=TRUE,
                                  BPPARAM=bpparam)
sampleInfo <- DataFrame(sample_name = sub(".bam", "", basename(bamFiles)),
                        file_bam = bamFiles,
                        group = rep(c("HinM", "Luc", "MinM"), each=6),
                        Treatment = factor(rep(rep(c("DOX", "NODOX"), each=3), 3),
                            levels=c("NODOX", "DOX")))
colData(C2C12.ens.SE) <- sampleInfo
save(C2C12.ens.SE, file=file.path(dataDir, "C2C12.ens.SE.rda"))

## DESeq2: HinM and MinM vs Luc
library(DESeq2)
se <- C2C12.ens.SE[rowSums(assay(C2C12.ens.SE)) > 9, ]
se.dox <- se[, colData(se)$Treatment == "DOX"]
C2C12.ens.ddsl <- lapply(c("HinM", "MinM"), function(x) {
    sub <- se.dox[, se.dox$group %in% c(x, "Luc")]
    print(dim(sub))
    sub$group <- factor(sub$group, levels=c("Luc", x))
    print(levels(sub$group))
    dds <- DESeqDataSet(sub, design = ~ group)
    dds <- DESeq(dds, parallel=TRUE)
})

names(C2C12.ens.ddsl) <- c("HinMvsLuc", "MinMvsLuc")
save(C2C12.ens.ddsl, file=file.path(pkgDir, "data", "C2C12.ens.ddsl.rda"))

## DESeq2: HinM and MinM vs NODOX
C2C12.ens.ddsl2 <- lapply(c("HinM", "MinM"), function(x) {
    sub <- se[, se$group==x]
    sub$Treatment <- factor(sub$Treatment, levels=c("NODOX", "DOX"))
    dds <- DESeqDataSet(sub, design = ~ Treatment)
    dds <- DESeq(dds, parallel=TRUE)
})
names(C2C12.ens.ddsl2) <- c("HinMvsNODOX", "MinMvsNODOX")
save(C2C12.ens.ddsl2, file=file.path(pkgDir, "data", "C2C12.ens.ddsl2.rda"))


#'
#' MB135 hg38 - MinM
#'
library(TxDb.Hsapiens.BioMart.ensembl.GRCh38.p2)
features <- exonsBy(TxDb.Hsapiens.BioMart.ensembl.GRCh38.p2, by="gene")
features <- keepStandardChromosomes(features, pruning.mode="fine")
seqlevelsStyle(features) <- "UCSC"

ngsDir <- "/shared/ngs/illumina/jwhiddon/151021_SN367_0569_AH7YMNBCXX"
bamDir <- file.path(ngsDir, "tophat", "bam")
bamFiles <- list.files(bamDir, pattern=".bam$", include.dirs=TRUE,
                       full.names=TRUE,
                       all.files=FALSE)
bamFiles <- bamFiles[1:6]
HinH.ens.SE <- summarizeOverlaps(features=features,
                                 reads=bamFiles,
                                 mode="IntersectionStrict",
                                 inter.feature=TRUE,
                                 singleEnd=TRUE,
                                 ignore.strand=TRUE,
                                 BPPARAM=bpparam)
sampleInfo <- DataFrame(sample_name = sub(".bam", "", basename(bamFiles)),
                         file_bam=bamFiles,
                         nick_name="HinH",
                         Treatment=factor(rep(c("NODOX", "DOX"), 3),
                             levels=c("NODOX", "DOX")))
colData(HinH.ens.SE) <- sampleInfo
save(HinH.ens.SE, file=file.path(dataDir, "HinH.ens.SE.rda"))

se <- HinH.ens.SE[rowSums(assay(HinH.ens.SE)) > 6, ]
HinH.ens.dds <- DESeqDataSet(se, design= ~ Treatment)
HinH.ens.dds <- DESeq(HinH.ens.dds, parallel=TRUE)
save(HinH.ens.dds, file=file.path(pkgDir, "data", "HinH.ens.dds.rda"))

