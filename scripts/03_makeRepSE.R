#!/app/R/devel/bin/Rscript
#./makeRepSE.R
#SBATCH -N1 -n6 -t 2-0 -p campus --mail-type=ALL --mail-user=cwon2@fhcrc.org -A tapscott_s

#'
#' This script consists with three components serving different purposes of
#' the analysis. It starts with hit-count for repNames followed by repeat enrichment
#' analyis (rmskStats) and finally, visualization.
#' 

library(GenomicAlignments)
library(canFam3.rmsk)
library(rmskStats)

#'
#' define parameters
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/canFam3.DuxFamily"
dataDir <- file.path(pkgDir, "data")
load(file.path(dataDir, "canine.ens.SE.rda"))
si <- colData(canine.ens.SE)[, c("sample_name", "file_bam", "pheno_type", "nick_name")]
workers <- cores <- 6L

#'
#' get repName counts
#' 
getRepNameFeatures <- function(rmsk) {
    features <- split(rmsk, rmsk$repName)

    #' fix (CATTC)n and (GAATC)n - split into two sattellite, simple_repeat
    tmp <- split(features[["(CATTC)n"]], factor(features[["(CATTC)n"]]$repClass))
    features[["(CATTC)n-Satellite"]] <- tmp[["Satellite"]]
    features[["(CATTC)n-Simple_repeat"]] <- tmp[["Simple_repeat"]]

    tmp <- split(features[["(GAATG)n"]], features[["(GAATG)n"]]$repClass)
    features[["(GAATG)n-Satellite"]] <- tmp[["Satellite"]]
    features[["(GAATG)n-Simple_repeat"]] <- tmp[["Simple_repeat"]]

    i <- which(names(features) %in% c("(CATTC)n", "(GAATG)n"))
    features <- features[-i]
    return(features)
}


canFam3.rmsk <- keepSeqlevels(canFam3.rmsk, value=seqlevels(canine.ens.SE),
                              pruning.mode="coarse")
features <- getRepNameFeatures(canFam3.rmsk)
all.reads <- rmskStats::sanitizeReads(bamFiles=si$file_bam, cores=workers)
save(all.reads, file=file.path(dataDir, "all.reads.rda"))

canine.repName.SE <- rmskStats::summarizeOverlaps.adjNH(features=features,
                                              all.reads=all.reads,
                                              type="within",
                                              ignore.strand=TRUE,
                                              inter.feature=TRUE)
colData(canine.repName.SE) <- si
rowData(canine.repName.SE) <-  as(t(sapply(features, rmskStats:::.getClassFamily)),
                            "DataFrame")
save(canine.repName.SE, file=file.path(dataDir, "canine.repName.SE.rda"))

cnt <- cbind(rowData(canine.repName.SE), as(assays(canine.repName.SE)[[1]], "DataFrame"))
write.csv(cnt, file.path(pkgDir, "inst", "stats", "rmsk_repName_Counts.csv"))
message("canine.repName.SE")

i <- cnt$repFamily == "Satellite"
cnt[i,]
#'
#' repeats analysis
#'                      
load(file.path(dataDir, "canine.repName.SE.rda"))
repStatsWrapper <- function(se, output.file, formula= ~ Treatment) {
    se <- se[rowSums(assays(se)[[1]]) > ncol(se), ]
    ddsRepName <- rmskStats:::repDESeq2(se, formula)
    sigRepName <- rmskStats:::repNameSig(ddsRepName, padj.thres=0.05, fc.thres=0.95)
    print(sigRepName)
    if (nrow(sigRepName) > 1) {
        repstats <- rmskStats:::repStats(ddsRepName, sigRepName)
        res <- rmskStats:::summarizeRepStats(repstats)
        res$total.sig.RepeatNames <- nrow(sigRepName)
        res$repstats <- repstats
        rmskStats:::repStatsReport(ddsRepName, sigRepName, repstats,
                               file=output.file, pval=0.05)
    } else { res <- "No significant repName."} 
    list(results=res, dds=ddsRepName)
}

rmsk.res <- lapply(c("HinC", "CinC", "CALTinC"), function(x) {
    idx <- canine.repName.SE$nick_name %in% c(x, "LinC")
    filename <- file.path(pkgDir, "inst", "stats", paste0("rmsk_", x, "_vs_LinC.xlsx"))
    se <- canine.repName.SE[, idx]
    se$nick_name <- factor(se$nick_name, levels=c("LinC", x))
    res <- repStatsWrapper(se, output.file=filename, formula= ~ nick_name)
})

rmsk.dds <- lapply(rmsk.res, "[[", "dds")
names(rmsk.dds) <- c("HinC", "CinC", "CALTinC")
save(rmsk.dds, file=file.path(dataDir, "rmsk.dds.rda"))

rmsk.res <- lapply(rmsk.res, "[[", "results")
names(rmsk.res) <- c("HinC", "CinC", "CALTinC")
save(rmsk.res, file=file.path(dataDir, "rmsk.res.rda"))

#'
#' print out the CALTinC statistics
#'
library(DESeq2)
calt <- rmsk.dds[["CALTinC"]]
res <- results(calt)
output <- counts(calt)
norm <- counts(calt, normalized=TRUE)
colnames(norm) <- paste0("norm_", colnames(norm))
output <- cbind(output, norm)
res <- append(mcols(rowRanges(calt))[, c("repClass", "repFamily")], res)
rownames(res) <- rownames(calt)
res <- append(res, as(output, "DataFrame"))
write.csv(as.data.frame(res),
          file=file.path(pkgDir, "inst", "stats", "rmsk_CALTinC_vs_LinC.csv"))

#'
#' visualization of enriched/depleted repeat families and classes
#' 
library(ggplot2)
library(DESeq2)
load(file.path(pkgDir, "data", "rmsk.dds.rda"))

#' (1) visuazlie significant repNames in terms of family and class
sig <- lapply(names(rmsk.dds)[c(1,2)], function(x, padj.thres=0.05, fc.thres=0.95) {
    ddsRepName <- rmsk.dds[[x]]
    res <- results(ddsRepName)
    res <- cbind(rowData(ddsRepName), res)
    keep <- which(res$padj < padj.thres & abs(res$log2FoldChange) > fc.thres)
    sig <- res[keep, ]
    sig$expression <- ifelse(sig$log2FoldChange > 0, "up", "down")
    sig$repFamily <- factor(sig$repFamily)
    sig$repClass <- factor(sig$repClass)
    sig$nick_name <- x
    sig <- as.data.frame(sig)
    title <- sprintf("%s: %3.0f and %3.0f repNames up- and down-regulated",
                     x,sum(sig$expression == "up"),
                     sum(sig$expression== "down"))
    ggfile <- file.path(figDir, paste0("rmsk_", x, "_SigRepName_on_repFamily.png"))
    gg <- ggplot(sig, aes(x=factor(repFamily), fill=expression)) +
    geom_bar(stat="count", width=0.7) +
        labs(title=title, x="repFamily", y="frequency of sig. repNames") +
            theme(axis.text.x=element_text(angle=90, hjust=1))
    ggsave(file=ggfile, gg)

    ggfile <- file.path(figDir, paste0("rmsk_", x, "_SigRepName_on_repClass.png"))
    gg <- ggplot(sig, aes(x=factor(repClass), fill=expression)) +
    geom_bar(stat="count", width=0.7) +
        labs(title=title, x="repClass", y="frequency of sig. repNames") +
            theme(axis.text.x=element_text(angle=90, hjust=1))
    ggsave(file=ggfile, gg)
    sig
})
names(sig) <- names(rmsk.dds)[c(1,2)]

#' (2) visualize enriched family and class
load(file.path(pkgDir, "data", "rmsk.res.rda"))
lapply(c("HinC", "CinC"), function(x) {
    repstats <- rmsk.res[[x]]$repstats
    ## family
    enriched_fam <- as.data.frame(repstats$fam_enrichment)
    enriched_fam <- subset(enriched_fam, num.sigRepName > 0)
    enriched_fam$repFamily <- rownames(enriched_fam)
    enriched_fam$enrichment <- ifelse(enriched_fam$prob < 0.05 &
                                          enriched_fam$num.sigRepName > enriched_fam$mu,
                                      "Yes", "No")
    ggfile <- file.path(figDir, paste0("rmsk_", x, "_enriched_repFam.png"))
    gg <- ggplot(enriched_fam, aes(x=repFamily, y=num.sigRepName, fill=enrichment)) +
        geom_bar(stat="identity", width=0.7) +
            geom_point(aes(x=repFamily, y=mu), color="red") + 
                theme(axis.text.x=element_text(angle=90, hjust=1)) +
                    labs(title=x, y="frequency of sig. repName") + 
                    scale_fill_manual(values=c("#999999", "#56B4E9")) 
    ggsave(file=ggfile, gg)

    ## class
    enriched_class <- as.data.frame(repstats$class_enrichment)
    enriched_class <- subset(enriched_class, num.sigRepName > 0)
    enriched_class$repClass <- rownames(enriched_class)
    enriched_class$enrichment <- ifelse(enriched_class$prob < 0.05 &
                                          enriched_class$num.sigRepName > enriched_class$mu,
                                      "Yes", "No")
    ggfile <- file.path(figDir, paste0("rmsk_", x, "_enriched_repClass.png"))
    gg <- ggplot(enriched_class, aes(x=repClass, y=num.sigRepName, fill=enrichment)) +
        geom_bar(stat="identity", width=0.7) +
            geom_point(aes(x=repClass, y=mu), color="red") + 
                theme(axis.text.x=element_text(angle=90, hjust=1)) +
                    labs(title=x, y="frequency of sig. repName") + 
                    scale_fill_manual(values=c("#999999", "#56B4E9")) 
    ggsave(file=ggfile, gg)
    invisible()
})

