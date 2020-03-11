# Function maps human seed ensembl ID to canine/mouse ensembl ID. Mapping could be
# 1:n. 

human_inparanoid_homology <- function(seed_ensembl) {
  require(DBI)
  require(hom.Hs.inp.db)
  require(dplyr)
  require(tidyr)
  require(purrr)
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  require(org.Cf.eg.db)

  human_canine <- as.list(hom.Hs.inpCANFA)
  human_mouse <- as.list(hom.Hs.inpMUSMU)

  # NOTE: seed_ensembl must be human ensembl; get protein IDs and symbol
   # (1) paired: human protein id -to- canine protein id (1:n)
  if (!is.character(seed_ensembl))
    stop("The seed ensembl ID must be a vector of character") 

  ensembl_prot <- 
    AnnotationDbi::select(org.Hs.eg.db, keys=seed_ensembl, 
                          columns=c("ENSEMBLPROT", "SYMBOL"),
                          keytype="ENSEMBL", multiVals="first") %>%
      dplyr::filter(!is.na(ENSEMBLPROT)) %>% #filter na that has no protein id
      dplyr::filter(ENSEMBLPROT %in% names(human_canine)) %>% # filter no homology in canine
      dplyr::filter(!duplicated(ENSEMBLPROT)) %>% # some protein ID have duplicated SYMBOL resulting duplicated protein ID
      dplyr::rename(HUMAN_ENSEMBL=ENSEMBL, HUMAN_SYMBOL=SYMBOL, HUMAN_ENSEMBLPROT=ENSEMBLPROT) # rename columns

  # (1) paired: human protein id -to- canine protein id (1:n)
  canine_paired <- mget(dplyr::pull(ensembl_prot, HUMAN_ENSEMBLPROT), hom.Hs.inpCANFA)   %>%
    tibble::enframe() %>% tidyr::unnest(cols=2) %>%
    dplyr::rename(HUMAN_ENSEMBLPROT=name, CANINE_ENSEMBLPROT=value)
  canine_prot <- # protein ID -to- ensembl ID/symbol; 1:n
    AnnotationDbi::select(org.Cf.eg.db, keys=canine_paired$CANINE_ENSEMBLPROT, 
                          columns=c("ENSEMBL", "SYMBOL"), # have duplicatd mapping
                          keytype="ENSEMBLPROT", multiVals="first") %>%
      dplyr::filter(!is.na(ENSEMBL)) %>% # remove no ensembl id mapping 
      dplyr::filter(!duplicated(ENSEMBL)) %>% # remove duplicated mapping 
      dplyr::rename(CANINE_ENSEMBL=ENSEMBL, CANINE_SYMBOL=SYMBOL, CANINE_ENSEMBLPROT=ENSEMBLPROT)
  # combind canine_paried and ensembl_prot => canine_prot  
  canine_paired <- canine_paired %>% dplyr::right_join(canine_prot, by="CANINE_ENSEMBLPROT")
  human_canine_paired <- ensembl_prot %>% right_join(canine_paired, by="HUMAN_ENSEMBLPROT")

  # (2) paired: human-to-mouse (1:n)
  ensembl_prot <- 
    AnnotationDbi::select(org.Hs.eg.db, keys=seed_ensembl, 
                          columns=c("ENSEMBLPROT", "SYMBOL"),
                          keytype="ENSEMBL", multiVals="first") %>%
      dplyr::filter(!is.na(ENSEMBLPROT)) %>% #filter na that has no protein id
      dplyr::filter(ENSEMBLPROT %in% names(human_mouse)) %>% # keep genes have homology in mouse
      dplyr::filter(!duplicated(ENSEMBLPROT)) %>% # some protein ID have duplicated SYMBOL resulting duplicated protein ID
      dplyr::rename(HUMAN_ENSEMBL=ENSEMBL, HUMAN_SYMBOL=SYMBOL, HUMAN_ENSEMBLPROT=ENSEMBLPROT) # rename columns

  mouse_paired <- mget(dplyr::pull(ensembl_prot, HUMAN_ENSEMBLPROT), hom.Hs.inpMUSMU)   %>%
    tibble::enframe() %>% tidyr::unnest(cols=2) %>%
    dplyr::rename(HUMAN_ENSEMBLPROT=name, MOUSE_ENSEMBLPROT=value)
  mouse_prot <- # protein ID -to- ensembl ID/symbol; 1:n
    AnnotationDbi::select(org.Mm.eg.db, keys=mouse_paired$MOUSE_ENSEMBLPROT, 
                          columns=c("ENSEMBL", "SYMBOL"), # have duplicatd mapping
                          keytype="ENSEMBLPROT", multiVals="first") %>%
      dplyr::filter(!is.na(ENSEMBL)) %>% # remove no ensembl id mapping 
      dplyr::filter(!duplicated(ENSEMBL)) %>% # remove duplicated mapping 
      dplyr::rename(MOUSE_ENSEMBL=ENSEMBL, MOUSE_SYMBOL=SYMBOL, MOUSE_ENSEMBLPROT=ENSEMBLPROT)    
  # combind canine_paried and ensembl_prot => canine_prot  
  mouse_paired <- mouse_paired %>% dplyr::right_join(mouse_prot, by="MOUSE_ENSEMBLPROT")
  human_mouse_paired <- ensembl_prot %>% right_join(mouse_paired, by="HUMAN_ENSEMBLPROT")
  
  return(list(human2canine=human_canine_paired, human2mouse=human_mouse_paired))
}