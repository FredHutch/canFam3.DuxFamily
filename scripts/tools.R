# lm_eqn is a copy from https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
lm_eqn <- function(df){
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)~"="~r2, 
         list(a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 2),
              r2 = format(sqrt(summary(m)$r.squared), digits = 3)))
    as.character(as.expression(eq));
}

z4state_hypergeometric <- function(z4state_homology_canine, res_DF) {
    ## hypergeometric testing for cleavage gene set (over-repersentation)
    res_df <- as(res_DF, "data.frame") %>%
      rownames_to_column(var="ENSEMBL")

    canine_z4state <- z4state_homology_canine %>%
      rename(ENSEMBL=CANINE_ENSEMBL) %>%
      left_join(res_df, by="ENSEMBL") %>%
      dplyr::filter(!is.na(log2FoldChange))

    universe <- pull(res_df, ENSEMBL)
    up_regulated <- res_df %>% dplyr::filter(log2FoldChange > 0 & padj < 0.05) %>%
      pull(ENSEMBL)
    black_ball <- length(universe) -  length(up_regulated) 
    white_ball <- length(up_regulated)
    ball_drown <- nrow(canine_z4state)
    white_drown <- canine_z4state %>%
      dplyr::filter(log2FoldChange > 0 & padj < 0.05) %>% nrow(.)
    prob <- dhyper(white_drown, m=white_ball, n=black_ball, k=ball_drown, log = FALSE)
    expected <- ball_drown * white_ball/(white_ball+black_ball)
    return(c(total=black_ball+white_ball, up_regulated=white_ball,
                  z4_set=ball_drown, up_reg_z4=white_drown,
                  expected=expected, prob=prob))
}

embryo_hypergeometric <- function(cleavage_homology_canine, res_DF) {
    ## hypergeometric testing for cleavage gene set (over-repersentation) or 
    ## human DUX4
    res_df <- as(res_DF, "data.frame") %>%
      rownames_to_column(var="ENSEMBL")
    canine_cleavage <- cleavage_homology_canine %>%
      rename(ENSEMBL=CANINE_ENSEMBL) %>%
      left_join(res_df, by="ENSEMBL") %>%
      dplyr::filter(!is.na(log2FoldChange))

    universe <- pull(res_df, ENSEMBL)
    up_regulated <- res_df %>% dplyr::filter(log2FoldChange > 0 & padj < 0.05) %>%
      pull(ENSEMBL)
    black_ball <- length(universe) -  length(up_regulated) 
    white_ball <- length(up_regulated)
    ball_drown <- nrow(canine_cleavage)
    white_drown <- canine_cleavage %>%
      dplyr::filter(log2FoldChange > 0 & padj < 0.05) %>% nrow(.)
    prob <- dhyper(x=white_drown, m=white_ball, n=black_ball, k=ball_drown, log = FALSE)
    expected <- ball_drown * white_ball/(white_ball+black_ball)
    return(c(total=black_ball+white_ball, up_regulated=white_ball,
                  gene_set=ball_drown, up_reg_geneset=white_drown,
                  expected=expected, prob=prob))
}

rmsk_enrichment <- function(rmsk_dds) {
    #dds <- rmsk.dds[[name]]
    rowdata <- as(rowData(rmsk_dds), "data.frame") %>%
      rownames_to_column(var="repName") %>%
      dplyr::select(repName, repClass, repFamily)
    # first threshold: p-value=0.0r correspondng to H_0: |lfc| < 0.5  
    res_df <- as(results(rmsk_dds, alpha=0.05, lfcThreshold=0.5), "data.frame") %>%
      rownames_to_column(var="repName") %>%
      left_join(rowdata, by="repName") 
    # padj: adj p-value =0.05 corresponding to H_0: |lfc| < 0.5
    # padj_null_zero: adj p-value=0.05 correspondng to H_0: lfc=0  
    universe <- nrow(res_df)
    selected_up <- res_df %>% 
      summarise(up = sum(padj < 0.05 & log2FoldChange > 0.5, na.rm=TRUE)) %>%
      pull(up)
    selected_down <- res_df %>% 
      summarise(down = sum(padj < 0.05 & log2FoldChange < 0.5, na.rm=TRUE)) %>%
      pull(down)

    # repFamily enrichment/depletion summary
    repFamily_summary <- res_df %>% group_by(repFamily) %>%
      summarise(total=n(), 
                up = sum(padj < 0.05 & log2FoldChange > 0, na.rm=TRUE), 
                down = sum(padj < 0.05 & log2FoldChange < 0, na.rm=TRUE),
                enrichment_prob = dhyper(x=up, k=total, m=selected_up, n=universe - selected_up),
                enrichment_mu = total * selected_up/universe,
                depletion_prob = dhyper(x=down, k=total, m=selected_down, n=universe-selected_down),
                depletion_mu = total * selected_down/universe)

    # repClass enrichement/depletion summary
    repClass_summary <- res_df %>% group_by(repClass) %>%
       summarise(total=n(), 
                up = sum(padj < 0.05 & log2FoldChange > 0, na.rm=TRUE), 
                down = sum(padj < 0.05 & log2FoldChange < 0, na.rm=TRUE),
                enrichment_prob = dhyper(x=up, k=total, m=selected_up, n=universe - selected_up),
                enrichment_mu = total * selected_up/universe,
                depletion_prob = dhyper(x=down, k=total, m=selected_down, n=universe-selected_down),
                depletion_mu = total * selected_down/universe)      
    # print out to a spread sheet                
    return(list(res_df=res_df, repFamily_summary=repFamily_summary, repClass_summary=repClass_summary))               
})

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