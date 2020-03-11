library(UCSCTrackTools)
pkgDir <- "~/tapscott/RNA-Seq/canFam3.DuxFamily"
load(file.path(pkgDir, "data", "canine.ens.SE.rda"))
df <- colData(canine.ens.SE)[, c("sample_name", "file_bam", "pheno_type", "nick_name")]
workers <- 4L

#'
#' coverage track
#' 
df$bwDir <- file.path("/fh/fast/tapscott_s/pub/tapscott/ucsc/bigWig",
                      paste0("canFam3_CnMb_", df$pheno_type))

library(RColorBrewer)
col <- col2rgb(brewer.pal(4, "Set2"))

lapply(c("LinC", "HinC", "CinC", "CALTinC"), function(x) {
    idx <- df$nick_name == x
    mycol <- col[, as.numeric(df$nick_name[idx])[1]]
    trackName <- paste0("CnMb_",df$pheno_type[idx][1], "_RNAseq")
    bwDir <- df$bwDir[idx][1]
    makeCoverageTracksWrapper(bamFiles=df$file_bam[idx],
                              bwDir=bwDir, #only one character
                              cores=3, trackName=trackName,
                              col=mycol)
})

#'
#' (1) make junction bed file
#' 
library(UCSCTrackTools)
library(GenomicAlignments)
pkgDir <- "~/tapscott/RNA-Seq/canFam3.DuxFamily"
load(file.path(pkgDir, "data", "canine.ens.SE.rda"))
df <- colData(canine.ens.SE)[, c("sample_name", "file_bam", "pheno_type", "nick_name")]
bedDir <- file.path(pkgDir, "inst", "junction_bed")
cores <- 4L
seqlev <- seqlevels(canine.ens.SE)
source("~/tapscott/R_package/UCSCTrackTools/R/makeJunctionBed.R")
source("~/tapscott/R_package/UCSCTrackTools/R/createHubTrackLine.bed.R")
bdFiles <- makeJunctionBed(bamFiles=df$file_bam, cores=cores, seqlev=seqlev,
                           min_junction_count=2,
                           genome="canFam3", verbose=TRUE,                   
                           ignore.strand=TRUE, outdir=bedDir) 
bdFiles <- list.files(bedDir, pattern="\\.bed$", full.names=TRUE)
names(bdFiles) <- sub(".bed", "", basename(bdFiles))
df$bedFiles <- bdFiles[rownames(df)] 
chrom_sizefile <- "/fh/fast/tapscott_s/CompBio/canFam3/canFam3.chrom.sizes"
group <- split(df, df$pheno_type)
#'
#' (2) and (3) convert bed to BigBed and create HubTrackLine
#' 
mclapply(group, function(x) {
    juncDir <- file.path("/fh/fast/tapscott_s/pub/tapscott/ucsc/junctions",
                         paste0("canFam3_CnMb_", x$pheno_type[1]))
    convertBedToBigBed(x$bedFiles, sample_name=NULL, cores=1L,
                       chrom_sizefile, outdir=juncDir)
    projectName <- basename(juncDir)
    trackName <- paste0(sub("canFam3_", "", projectName), "_junc")
    createHubTrackLine.junc.bed(bbDir=juncDir,
                                trackName=trackName,
                                pattern="\\.bb$")
}, mc.cores=cores, mc.preschedule=FALSE)


#'
#' testing makeJunctionWrapper()
#'
pkgDir <- "/fh/fast/tapscott_s/R/RNA-Seq/canFam3.DuxFamily"
load(file.path(pkgDir, "data", "canine.ens.SE.rda"))
ucscDir <- "/fh/fast/tapscott_s/pub/tapscott/ucsc/junctions"
ngsDir <- "/shared/ngs/illumina/atyler/170208_SN367_0857_AHF37NBCXY"
bamDir <- file.path(ngsDir, "tophat", "bam")
bamFiles <- list.files(bamDir, pattern="\.bam$", full.names=TRUE)[c(1,5, 6)]
seqlev <- seqlevels(canine.ens.SE)
chrom_sizefile <- "/fh/fast/tapscott_s/CompBio/canFam3/canFam3.chrom.sizes"
bedDir <- file.path(pkgDir, "inst", "junction_bed")
bigBedDir <- file.path(ucscDir, "canFam3_CnMb_luciferase")
trackName <- "CnMb_luciferase_junc"

makeJunctionTracksWrapper(bamFiles, genome="canFam3",
                          cores=1L, seqlev=seqlev,
                          min_junction_count=2,
                          chrom_sizefile=chrom_sizefile,
                          bedDir=bedDir,
                          bigBedDir=bigBedDir,
                          trackName=trackName)
                          
                          
    
