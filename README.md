## Genome build
- UCSC Genome broswser assembly ID: canFam3
- NCBI Genome information: NCBI genome/85 CanFam3.1 GCF_000002285.3

software: Bioconductor 3.5 / R 3.4

## Annotation
- Bioconductor TxDb package: Do not use _TxDb.Cfamiliaris.UCSC.canFam3.refGene_. The refGene track on UCSC is not correct. 
- Ensembl86 for canFam3.1 is also available on UCSC. We have built an Ensembl 
  TxDb package:/fh/fast/tapscott_s/CompBio/canFam3/TxDb.Cfamiliaris.UCSC.canFam3.ensGene_1.0.0.tar.gz.
- Bioconductor BSgenome package: _BSgenome.Cfamiliaris.UCSC.canFam3_
- Bioconductor Genome wide annotation package: _org.Cf.eg.db_
- Homology: downloaded from Ensembl and save at `~/CompBio/R_package/ensOrtholog.db/data/mouse_human_ensOrtholog.rda`
- RMSK: download from UCSC genome browser canFam3 and make
_canFam3.rmsk_ pacakge.  The package is stored at
`~/CompBio/canFam3/canFam3_1.0.0.tar.gz`.

## Where is row data
/shared/ngs/illumina/atyler/170208_SN367_0857_AHF37NBCXY/Unaligned/

## Standard chromosomes
Use standard chromosomes
seqlevels of canFam3:
 sl[!grepl("chrUn_", sl)]
 [1] "chr1"  "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17"
[10] "chr18" "chr19" "chr2"  "chr20" "chr21" "chr22" "chr23" "chr24" "chr25"
[19] "chr26" "chr27" "chr28" "chr29" "chr3"  "chr30" "chr31" "chr32" "chr33"
[28] "chr34" "chr35" "chr36" "chr37" "chr38" "chr4"  "chr5"  "chr6"  "chr7"
[37] "chr8"  "chr9"  "chrM"  "chrX"

## Differential Analysis
**step 1**
features: TxDb.Cfamiliaris.UCSC.canFam3.ensGene
count: summarizedExperiment
differential analysis: CinC vs Luc, HinC vs Luc, CALTinC
hypothesis: mean log fold change greater than 1.
repeat analysis: need to be re-done by the countRepeat package.

