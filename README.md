![duxc-summary]("./gitbook/images/duxc_summary/png")

To support the reproducibility and transparency of the computational work for the manuscript, [_Canine DUXC: Implications for DUX4 retrotransposition and preclinical models of FSHD_](https://academic.oup.com/hmg/advance-article/doi/10.1093/hmg/ddab352/6457948), we made this repository to include our processed RNA-seq and ChIP-seq datasets, as well as a [__book__](https://fredhutch.github.io/canFam3.DuxFamily/) containing nine chapters demonstrating all the bioinformatic analysis and automatically executable R code (you can directly run the code and reproduce the statistics and figures in the manuscript).


## Datasets 
Datasets are stored in the __/data__ folder.

- canine RNA-seq ([GSE188928](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188928)): DUX4, DUXC, and DUX-ALT transcriptome in canine skeletal muscle and repeat elements expression stored in `DESeqDataSet` instances
```
├── CALTinC.ens.dds.rda: DUXC-ALT transcriptome in canine skeletal muscle
├── CinC.ens.dds.rda: DUXC transcriptome
├── HinC.ens.dds.rda: DUX4 transcriptome in canine skeletal muscle
├── rmsk.dds.rda: DUX4, DUXC, DUXC-ALT repeat elements expression
├── CinC_rmsk.rda: DUXC repeat elements DESeq2 results
├── HinC_rmsk.rda: DUX4 repeat elements DESeq2 results
├── CALTinC_rmsk.rda: DUXC-ALT repeat elements DESeq2 results
```

- human DUX4 transcriptome in human myoblast ([GSE85461](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85461)) and mouse DUX transcriptome in mouse myoblast cell lines ([GSE87282](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87282))     
```
├── C2C12.ens.ddsl2.rda: DUX and DUX4 transcriptome in mouse myoblast cell lines
├── HinH.ens.dds.rda: DUX4 transcriptome in human myoblast
├── HinH_rmsk.rda: repeat element expression in HinH transcriptome (DESeq2 results)
```

- miscellaneous. The orthologue are collected from Ensembl orthologue database V88.   
```
├── human_cleavage.rda: Hendrikson 2017 zygotic activation gene set
├── canine_human_mouse_ensOrtholog.rda: human-canine and mouse-canine paired orthologues
├── human_mouse_ensOrtholog.rda: human-mouse paired Orthologue
├── peaks_list.rda: ChIP-seq peaks 
├── chipseq_si.rda: ChIP-seq sample information
└── z4_ensembl.rda: Akiyama 2015 2C-like gene set from 
```

## Genome build / software / annotation

- UCSC Genome browser assembly ID: canFam3       
- NCBI Genome information: NCBI genome/85 CanFam3.1 GCF_000002285.3     

__Software__: Bioconductor 3.10 / R 3.6.3

__Annotation__:

- Ensembl86 for canFam3.1 is also available on UCSC. We have built a Bioconductor TxDb package:`/fh/fast/tapscott_s/CompBio/canFam3/TxDb.Cfamiliaris.UCSC.canFam3.ensGene_1.0.0.tar.gz`. This [page](https://fredhutch.github.io/canFam3.DuxFamily/) gives codes for acquiring Ensembl database and building a TxDb package.   
- Bioconductor BSgenome package: _BSgenome.Cfamiliaris.UCSC.canFam3_     
- Bioconductor Genome wide annotation package: _org.Cf.eg.db_     
- Homology: downloaded from Ensembl and saved as    `/mouse_human_ensOrtholog.rda`
- RMSK: download from UCSC genome browser canFam3 and made an R
_canFam3.rmsk_ package. 

## Where is the row data

- GEO superseries (public): [GSE188928](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188928)
- FredHutch share drive (private): `/shared/ngs/illumina/atyler/170208_SN367_0857_AHF37NBCXY/Unaligned/`





