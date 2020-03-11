#!/app/R/devel/bin/Rscript
#./makeSE.R
#SBATCH -N1 -n12 -t 2-0 -p campus --mail-type=ALL --mail-user=cwon2@fhcrc.org -A tapscott_s

pkgDir <- "~/tapscott/RNA-Seq/canFam3.DuxFamily"
dataDir <- file.path(pkgDir, "data")
ngsDir <- "/shared/ngs/illumina/atyler/170208_SN367_0857_AHF37NBCXY"
bamDir <- file.path(ngsDir, "tophat", "bam")
bamFiles <- list.files(bamDir, pattern="\\.bam$", full.names=TRUE)
bamFiles <- bamFiles[c(1, 5, 6, 7, 8, 9, 10, 11, 12, 2, 3, 4)]

#'
#' libraries and functions
#' 
library(GenomicFeatures)
library(GenomicAlignments)
library(TxDb.Cfamiliaris.UCSC.canFam3.refGene)
library(TxDb.Cfamiliaris.UCSC.canFam3.ensGene)
library(org.Cf.eg.db)
source("~/tapscott/R_package/SeqPipelineTools/R/makeSE.R")
source("~/tapscott/R_package/SeqPipelineTools/R/simpleDESeqReport.R")
workers <- 12L

#'
#' feature prediction
#'
ref_exons <- exonsBy(TxDb.Cfamiliaris.UCSC.canFam3.refGene, by="gene")
ref_exons <- keepStandardChromosomes(ref_exons, pruning.mode="fine")
ens_exons <- exonsBy(TxDb.Cfamiliaris.UCSC.canFam3.ensGene, by="gene")
ens_exons <- keepStandardChromosomes(ens_exons, pruning.mode="fine")

#'
#' make SummearizedExperiments for refGene and ensGene
#'
nick_name <- c(rep("LinC", 3), rep("HinC", 3),
               rep("CinC", 3), rep("CALTinC", 3))
nick_name <- factor(nick_name, levels=c("LinC", "HinC", "CinC", "CALTinC"))
pheno_type <- sapply(strsplit(basename(bamFiles), "_"), "[[", 2)
pheno_type <- factor(pheno_type, levels=c("luciferase", "hDUX4CA", "cDUXCintac", "cDUXCinter"))

#' refGene track
canine.ref.SE <- makeSEwrapper(bamFiles, pkgDir=pkgDir, OrgDb=org.Cf.eg.db,
                           TxDb=NULL, features=ref_exons, workers=workers,
                           pheno_type=nick_name, title="Canine_refGene",
                           plot.PCA=TRUE, gene_id_type="ENTREZID",
                           do.DESeq=FALSE, print.report=FALSE)$se
canine.ref.SE$nick_name  <- canine.ref.SE$pheno_type
canine.ref.SE$pheno_type <-  pheno_type
save(canine.ref.SE, file=file.path(dataDir, "canine.ref.SE.rda"))

#' ensGene track
canine.ens.SE <- makeSEwrapper(bamFiles, pkgDir=pkgDir, OrgDb=org.Cf.eg.db,
                           TxDb=NULL, features=ens_exons, workers=workers,
                           pheno_type=nick_name, title="Canine_ensGene",
                           plot.PCA=TRUE, 
                           do.DESeq=FALSE, print.report=FALSE)$se
canine.ens.SE$nick_name  <- canine.ens.SE$pheno_type
canine.ens.SE$pheno_type <-  pheno_type
save(canine.ens.SE, file=file.path(dataDir, "canine.ens.SE.rda"))

#'
#' refGene: three different DESeq analysis
#' 
sel <- lapply(c("HinC", "CinC", "CALTinC"), function(x) {
    idx <- nick_name %in% c(x, "LinC")
    title <- paste0("ref_", x, "_vs_LinC")
    factor <- factor(nick_name[idx], levels=c("LinC", x))
    res <- makeSEwrapper(bamFiles[idx], pkgDir=pkgDir, OrgDb=org.Cf.eg.db,
                        TxDb=NULL, features=ref_exons, workers=workers,
                        pheno_type=factor, title=title,
                        plot.PCA=FALSE, lfcThreshold=1,
                        do.DESeq=TRUE, print.report=TRUE)
    se <- res$se
    se$nick_name  <- se$pheno_type
    se$pheno_type <-  pheno_type[idx]
    res$se <- se
    return(res)
})
names(sel) <- c("HinC", "CinC", "CALTinC")
HinC.ref.dds    <- sel[["HinC"]]$dds
CinC.ref.dds    <- sel[["CinC"]]$dds
CALTinC.ref.dds <- sel[["CALTinC"]]$dds
save(HinC.ref.dds,    file=file.path(dataDir, "HinC.ref.dds.rda"))
save(CinC.ref.dds,    file=file.path(dataDir, "CinC.ref.dds.rda"))
save(CALTinC.ref.dds, file=file.path(dataDir, "CALTinC.ref.dds.rda"))

#'
#' ensGene: three different DESeq analysis
#' 
sel <- lapply(c("HinC", "CinC", "CALTinC"), function(x) {
    idx <- nick_name %in% c(x, "LinC")
    title <- paste0("ens_", x, "_vs_LinC")
    factor <- factor(nick_name[idx], levels=c("LinC", x))
    res <- makeSEwrapper(bamFiles[idx], pkgDir=pkgDir, OrgDb=org.Cf.eg.db,
                         TxDb=NULL, features=ens_exons, workers=workers,
                         pheno_type=factor, title=title,
                         plot.PCA=FALSE, lfcThreshold=1,
                         gene_id_type="ENSEMBL",
                         do.DESeq=TRUE, print.report=TRUE)
    se <- res$se
    se$nick_name  <- se$pheno_type
    se$pheno_type <-  pheno_type[idx]
    res$se <- se
    return(res)
})
names(sel) <- c("HinC", "CinC", "CALTinC")
HinC.ens.dds    <- sel[["HinC"]]$dds
CinC.ens.dds    <- sel[["CinC"]]$dds
CALTinC.ens.dds <- sel[["CALTinC"]]$dds
save(HinC.ens.dds,    file=file.path(dataDir, "HinC.ens.dds.rda"))
save(CinC.ens.dds,    file=file.path(dataDir, "CinC.ens.dds.rda"))
save(CALTinC.ens.dds, file=file.path(dataDir, "CALTinC.ens.dds.rda"))

