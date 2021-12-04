To support the reproducibility and transparancy of the computational work for the manuscript, _Canine DUXC: Implications for DUX4 retrotransposition and preclinical models of FSHD_, we made this repository to include our processed RNA-seq and ChIP-seq dataset, as well as our preprocess pipeline and analysis.  We also made a gitbook containning nine chapters demonstrating the analysis workflow and automatically executable R code (you can directly run the code and reproduce the statistics and figures in the manuscript).


## Genome build / software / annotation

- UCSC Genome broswser assembly ID: canFam3       
- NCBI Genome information: NCBI genome/85 CanFam3.1 GCF_000002285.3     

__Software__: Bioconductor 3.5 / R 3.4

__Annotation__:

- Bioconductor TxDb package: Do not use _TxDb.Cfamiliaris.UCSC.canFam3.refGene_. The refGene track on UCSC is not correct. Use standard chromosomes.   
- Ensembl86 for canFam3.1 is also available on UCSC. We have built an Ensembl 
  TxDb package:/fh/fast/tapscott_s/CompBio/canFam3/TxDb.Cfamiliaris.UCSC.canFam3.ensGene_1.0.0.tar.gz.     
- Bioconductor BSgenome package: _BSgenome.Cfamiliaris.UCSC.canFam3_     
- Bioconductor Genome wide annotation package: _org.Cf.eg.db_     
- Homology: downloaded from Ensembl and save at      `~/CompBio/R_package/ensOrtholog.db/data/mouse_human_ensOrtholog.rda`
- RMSK: download from UCSC genome browser canFam3 and make
_canFam3.rmsk_ pacakge.  The package is stored at
`~/CompBio/canFam3/canFam3_1.0.0.tar.gz`.     

## Where is row data
FredHutch share drive: `/shared/ngs/illumina/atyler/170208_SN367_0857_AHF37NBCXY/Unaligned/`





