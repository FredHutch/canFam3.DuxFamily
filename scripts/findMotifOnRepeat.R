library(rtracklayer)
motifGff <- "~/tapscott/canFam3/HinM_MEME_FIMO_locationsIn_CanFam3.gff"
h_motif <- import(motifGff)
library(canFam3.rmsk)
library(rmskStats)
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/canFam3.DuxFamily"
ov <- findOverlaps(canFam3.rmsk, h_motif, minoverlap=10L, ignore.strand=TRUE)
canFam3.rmsk$motif_hit <- FALSE
canFam3.rmsk$motif_id <- NA
canFam3.rmsk$motif_hit[queryHits(ov)] <- TRUE
canFam3.rmsk$motif_id[queryHits(ov)] <- h_motif$ID[subjectHits(ov)]

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
    mcols(features) <- t(sapply(features, function(x) {
        c(repFamily=as.character(x$repFamily)[1],
          repClass=as.character(x$repClass)[1])
    }))
    
    return(features)
}

rmsk <- getRepNameFeatures(canFam3.rmsk)

df <- DataFrame(repName=names(rmsk),
                num.location=elementNROWS(rmsk),
                num.motif_hits=sapply(rmsk, function(x) sum(x$motif_hit)),
                repFamily=mcols(rmsk)$repFamily,
                repClass=mcols(rmsk)$repClass)
               
df$hit_freq <- df$num.motif_hits/df$num.location
df <- df[order(df$hit_freq, decreasing=TRUE), ]
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/canFam3.DuxFamily"
write.csv(df, file=file.path(pkgDir, "inst", "stats", "RepName_human_motif_hits.csv"))

## why (ATCG)n - simple repeat, (CGATG)n- simple repeat, HSATII -satellite
## have no repClass and repFamily

tmp <- df[df$num.motif_hits > 0, ]
