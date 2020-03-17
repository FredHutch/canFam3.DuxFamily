# lm_eqn is a copy from https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
lm_eqn <- function(df){
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)~"="~r2, 
         list(a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 2),
              r2 = format(sqrt(summary(m)$r.squared), digits = 3)))
    as.character(as.expression(eq));
}

gsea_hypergeomatric <- function(universe, de_genes, gene_set) {
    # this function performs hypergeometric testing for Gene Set Enrichment Analysis
}

.do_goseq <- function(universe, selected, threshold_pval=0.01, 
                      return.DEInCat=FALSE, dds=NULL) {
    # This function perfroms GO term analysis for canFam3 genome build genes                          
    require(goseq)
    require(org.Cf.eg.db)
    require(GO.db)
    require(TxDb.Cfamiliaris.UCSC.canFam3.ensGene)

    txsByGene <- transcriptsBy(TxDb.Cfamiliaris.UCSC.canFam3.ensGene, by="gene")
    names(txsByGene) <- 
        sapply(strsplit(names(txsByGene), ".", fixed=TRUE), "[[", 1)
    lengthData <- median(width(txsByGene))
    
    isDEGs <- as.integer(universe %in% selected)
    names(isDEGs) <- universe
    bias.data <- lengthData[names(isDEGs)]
    pwf <- nullp(isDEGs, bias.data=bias.data, plot.fit=FALSE)
    # rtracklayer::ucscGenomes()[40:42] does support canFam3
    # but supportedGeneIDs() does onot support canFam3??
    GO.BP <- goseq(pwf, "canFam3", "ensGene", test.cats=c("GO:BP"))
    
    enriched.BP <- GO.BP %>%
      mutate(fdr = p.adjust(over_represented_pvalue, method="BH")) %>%
      dplyr::filter(over_represented_pvalue < threshold_pval)
    
    if (return.DEInCat & !is.null(dds)) {
      cat_genes <- lapply(enriched.BP$category, function(GOID) {
        cat_genes <- .cat2DEgenes(GOID, pwf=pwf)
        cat_genename <- .mapID2GeneName(dds, 
                                        id=cat_genes, clean=TRUE)
        paste(cat_genename, collapse=",")
      })  
      enriched.BP <- add_column(enriched.BP, DEInCat=cat_genes)
    }

    return(enriched.BP)
}

.cat2DEgenes <- function(GOID, pwf) {
    gene2cat <- getgo(rownames(pwf), "hg38", "ensGene", fetch.cats = "GO:BP")
    names(gene2cat) <- rownames(pwf)
    cat2gene <- goseq:::reversemapping(gene2cat)
    #' sanity check
    doesIDexist <- GOID %in% names(cat2gene)
    if (!all(doesIDexist)) stop("GOID is not found")
    sig_gene <- rownames(pwf)[pwf$DEgenes==1]
    geneInGOID <- cat2gene[[GOID]]
    sig_gene[sig_gene %in% geneInGOID]
}

.mapID2GeneName <- function(dds, id, clean=TRUE) {
    if (clean)  {
      rownames(dds) <- sapply(strsplit(rownames(dds), ".", fixed=TRUE),
                    "[[", 1)
    }
    rowData(dds[id])$gene_name
}