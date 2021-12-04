To support the reproducibility and transparancy of the computational work for the manuscript, _Canine DUXC: Implications for DUX4 retrotransposition and preclinical models of FSHD_, we made this repository to include our processed RNA-seq and ChIP-seq dataset, as well as our preprocess pipeline and analysis.  We also made a gitbook containning nine chapters demonstrating the analysis workflow and automatically executable R code (you can directly run the code and reproduce the statistics and figures in the manuscript).


## Datasets

- canine RNA-seq [GSE188928](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188928): DUX4, DUXC, and DUX-ALT transcriptome in canine skeletal muscle and repeat elements expression stored in `DESeqDataSet` instances
```
├── CALTinC.ens.dds.rda: DUXC-ALT transcriptome in canine skeletal muscle
├── CinC.ens.dds.rda: DUXC transcriptome
├── HinC.ens.dds.rda: DUX4 transcriptome in canine skeletal muscle
├── rmsk.dds.rda: DUX4, DUXC, DUXC-ALT repeat elements expression
├── CinC_rmsk.rda: DUXC repeat elements DESeq2 results
├── HinC_rmsk.rda: DUX4 repeat elements DESeq2 resutls
├── CALTinC_rmsk.rda: DUXC-ALT repeat elements DESeq2 results
```

- human DUX4 transcriptome in human myblast([GSE85461](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85461)) and mouse Dux transcriptome in mouse myblast cell lines ([GSE87282](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87282))     
```
├── C2C12.ens.ddsl2.rda: Dux and DUX4 transcriptome in mouse myoblast cell lines
├── HinH.ens.dds.rda: DUX4 transcriptome in human myoblast
├── HinH_rmsk.rda: repeat element expression in HinH transcriptome (DESeq2 resutls)
```

- miscellaneous. The orthologous are collected from Ensembl ortholog database V88.   
```
├── human_cleavage.rda: Hendrikson 2017 zygotic activation gene set
├── canine_human_mouse_ensOrtholog.rda: human-canine and mouse-canine paired orthologous
├── human_mouse_ensOrtholog.rda: human-mouse paired Ortholog
├── peaks_list.rda: ChIP-seq peaks 
├── chipseq_si.rda: ChIP-seq sample information
└── z4_ensembl.rda: Akiyama 2015 2C-like gene set from 
```

## Genome build / software / annotation

- UCSC Genome broswser assembly ID: canFam3       
- NCBI Genome information: NCBI genome/85 CanFam3.1 GCF_000002285.3     

__Software__: Bioconductor 3.10 / R 3.6.3

__Annotation__:

- Ensembl86 for canFam3.1 is also available on UCSC. We have built a Bioconductor TxDb package:`/fh/fast/tapscott_s/CompBio/canFam3/TxDb.Cfamiliaris.UCSC.canFam3.ensGene_1.0.0.tar.gz`. This [page](https://fredhutch.github.io/canFam3.DuxFamily/) gives codes for acquiring Ensembl database and building a TxDb package.   
- Bioconductor BSgenome package: _BSgenome.Cfamiliaris.UCSC.canFam3_     
- Bioconductor Genome wide annotation package: _org.Cf.eg.db_     
- Homology: downloaded from Ensembl and saved as    `/mouse_human_ensOrtholog.rda`
- RMSK: download from UCSC genome browser canFam3 and made an R
_canFam3.rmsk_ pacakge. 

## Where is the row data

- GEO superseries (public): [GSE188928](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188928)
- FredHutch share drive (private): `/shared/ngs/illumina/atyler/170208_SN367_0857_AHF37NBCXY/Unaligned/`





