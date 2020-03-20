# this script is to study (1) DUXC transicription activation, (2) GSEA enrichement
# (3) DUC4 trasncription activation in canine myoblast and (4) repeat family enrichment 
# anlaysis


pkg_dir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/canFam3.DuxFamily"
fig_dir <- file.path(pkg_dir, "figures")
data_dir <- file.path(pkg_dir, "data")
source(file.path(pkg_dir, "scripts", "inparanoid_homology.R"))
source(file.path(pkg_dir, "scripts", "tools.R"))
# loading library
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(org.Cf.eg.db)

# loading datasets
load(file.path(data_dir, "human_cleavage.rda"))
load(file.path(data_dir, "CinC.ens.dds.rda"))
load(file.path(data_dir, "CALTinC.ens.dds.rda"))
load(file.path(data_dir, "HinC.ens.dds.rda"))
load(file.path(data_dir, "HinH.ens.dds.rda"))
rownames(HinC.ens.dds) <- sapply(strsplit(rownames(HinC.ens.dds), fixed=TRUE, "."), "[", 1)
rownames(CinC.ens.dds) <- sapply(strsplit(rownames(CinC.ens.dds), fixed=TRUE, "."), "[", 1)
rownames(CALTinC.ens.dds) <- sapply(strsplit(rownames(CALTinC.ens.dds), fixed=TRUE, "."), "[", 1)

#
# (0) tidy datasets: DESeq2 results
#

# DataFrame form
CinC_res <- results(CinC.ens.dds, lfcThreshold=1, alpha=0.05)
HinC_res <- results(HinC.ens.dds, lfcThreshold=1, alpha=0.05)
CALTinC_res <- results(CALTinC.ens.dds, lfcThreshold=1, alpha=0.05)

# data.frame form
CinC_res.2 <- as(CinC_res, "data.frame") %>%
  rownames_to_column(var="ENSEMBL") %>%
  dplyr::select(ENSEMBL, log2FoldChange, padj) %>%
  rename(CinC_logFC=log2FoldChange, CinC_padj=padj) %>%
  dplyr::mutate(CinC_DE=CinC_padj < 0.05)

HinC_res.2 <- as(HinC_res, "data.frame") %>%
  rownames_to_column(var="ENSEMBL") %>%
  dplyr::select(ENSEMBL, log2FoldChange, padj) %>%
  rename(HinC_logFC=log2FoldChange, HinC_padj=padj) %>%
  dplyr::mutate(HinC_DE=HinC_padj < 0.05)

CALTinC_res.2 <- as(CALTinC_res, "data.frame") %>%  
  rownames_to_column(var="ENSEMBL") %>%
  dplyr::select(ENSEMBL, log2FoldChange, padj) %>%
  rename(CALTinC_logFC=log2FoldChange, CALTinC_padj=padj) %>%
  dplyr::mutate(CALTinC_DE=CALTinC_padj < 0.05)


#####################################
#
# (1) DESeq2 resutls and MAplots
#
####################################
# (a) CinC
summary(CinC_res)
pdf(file.path(fig_dir, "CinC_maplot.pdf"), width=4, height=4)
plotMA(CinC_res, alpha=0.05, main="CinC")
text(x=10000, y=5, labels="up: 1562", adj=c(0, 0))
text(x=10000, y=-5, labels="down: 191", adj=c(0, 1))
dev.off()

#' how many cleavege genes are up-regulated?
#' do down-regulated genes overlap between three senarios?

# (b) HinC
summary(HinC_res)
pdf(file.path(fig_dir, "HinC_maplot.pdf"), width=4, height=4)
plotMA(HinC_res, alpha=0.05, main="HinC")
text(x=10000, y=5, labels="up: 1553", adj=c(0, 0))
text(x=10000, y=-5, labels="down: 231", adj=c(0, 1))
dev.off()

# (c) CALTinC
summary(CALTinC_res)
pdf(file.path(fig_dir, "CALTinC_maplot.pdf"), width=4, height=4)
plotMA(CALTinC_res, alpha=0.05, main="CALTinC")
text(x=10000, y=2, labels="up: 20", adj=c(0, 0))
text(x=10000, y=-2, labels="down: 114", adj=c(0, 1))
dev.off()

###################################
#
# (2) comparison by scatter plot?
#
##################################
# (a) CinC vs CALT
CALTinC_CinC <- CinC_res.2 %>% inner_join(CALTinC_res.2, by="ENSEMBL") %>%
  dplyr::mutate(both_DE= CALTinC_DE & CinC_DE) %>%
  dplyr::mutate(DE_status = "neither") %>%
  dplyr::mutate(DE_status = replace(DE_status, CinC_DE, "either")) %>%
  dplyr::mutate(DE_status = replace(DE_status, CALTinC_DE, "either")) %>%
  dplyr::mutate(DE_status = replace(DE_status, both_DE, "both")) %>%
  dplyr::mutate(DE_status = factor(DE_status, levels=c("neither", "either", "both")))

table(CALTinC_CinC$DE_status)

gg <- ggplot(CALTinC_CinC, aes(x=CinC_logFC, y=CALTinC_logFC)) +
  geom_point(size=0.5, alpha=0.3, aes(color=DE_status)) +
  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
  labs(title="CinC vs CALTinC logFC", color="DE status") +
  theme_bw()+
  theme(legend.position = c(0.13, 0.825) )
pdf(file.path(pkg_dir, "figures", "CALTinC_CinC_logFC_scatter.pdf"),
    height=4, width=4)
plot(gg)  
dev.off()    

# (b) CinC vs HinC
HinC_CinC <- HinC_res %>%
  left_join(CinC_res, by="ENSEMBL") %>%
  dplyr::mutate(both_DE= HinC_DE & CinC_DE) %>%
  dplyr::filter(!is.na(HinC_logFC), !is.na(CinC_logFC)) %>%
  dplyr::mutate(DE_status = "neither") %>%
  dplyr::mutate(DE_status = replace(DE_status, HinC_DE, "either")) %>%
  dplyr::mutate(DE_status = replace(DE_status, CinC_DE, "either")) %>%
  dplyr::mutate(DE_status = replace(DE_status, both_DE, "both")) %>%
  dplyr::mutate(DE_status = factor(DE_status, levels=c("neither", "either", "both")))

summary <-  HinC_CinC %>%
  summarise(both=sum(both_DE, na.rm=TRUE),
            HinC=sum(HinC_DE, na.rm=TRUE), 
            CinC=sum(CinC_DE, na.rm=TRUE),
            cor = cor(HinC_logFC, CinC_logFC))

gg <- ggplot(HinC_CinC, aes(x=CinC_logFC, y=HinC_logFC)) +
  geom_point(size=0.5, alpha=0.3, aes(color=DE_status)) +
  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
  geom_smooth(method="lm", se=FALSE, show.legend=FALSE, color="gray70", alpha=0.5,
              linetype="dashed") +
  labs(title="HinC vs CinC logFC", color="DE status") +
  theme_bw()+
  theme(legend.position = c(0.15, 0.8)) 
df <- HinC_CinC %>% dplyr::select(HinC_logFC, CinC_logFC) %>%
  rename(y=HinC_logFC, x=CinC_logFC)  
gg <- gg + geom_text(x = 0, y = -5, vjust=1, hjust=0, 
  label = lm_eqn(df), parse = TRUE, color="gray50") 
pdf(file.path(pkg_dir, "figures", "HinC_CinC_logFC_scatter.pdf"),
    height=4, width=4)
plot(gg)  
dev.off()  

######################################
#
# (3) cleavage-stage gene signature
#
######################################
cleavage_homology <- human_inparanoid_homology(as.character(human_cleavage$GeneID))
canine_cleavage <- cleavage_homology$human2canine %>%
  rename(ENSEMBL=CANINE_ENSEMBL) %>%
  left_join(CinC_res.2, by="ENSEMBL") %>%
  left_join(HinC_res.2, by="ENSEMBL") # total 73 genes
message("Total", nrow(canine_cleanve), "human-to-canine homologues.")

# use functions from tools.R
CinC_enrichment_prob <- 
  geneset_hypergeometric(cleavage_homology$human2canine, CinC_res)
HinC_enrichement_prob <- 
  geneset_hypergeometric(cleavage_homology$human2canine, HinC_res)
knitr::kable(rbind(CinC_enrichment_prob, HinC_enrichement_prob))  


######################################
#
# (4) mouse 2C-state gene set enrichment
#
######################################
# test: mouse to caninef
library(xlsx)
library(org.Mm.eg.db)
source(file.path(pkg_dir, "scripts", "inparanoid_homology.R"))

z4 <- read.xlsx(file.path(pkg_dir, "extdata", "ZSCAN4_ES_Upreg_Fixed.xlsx"),
                      sheetIndex=1, stringsAsFactors=FALSE)
z4_ensembl <- AnnotationDbi::select(org.Mm.eg.db, 
                             keys=as.character(z4$GeneID), columns=c("ENSEMBL"), 
                             keytype="REFSEQ") 
                               
                       
seed_ensembl <- unique(z4_ensembl$ENSEMBL)
mouse2canine <- mouse_inparanoid_homology(seed_ensembl)
dim(mouse2canine)
CinC_enrichment_prob <- 
  z4state_hypergeometric(mouse2canine, CinC_res)
HinC_enrichement_prob <- 
  z4state_hypergeometric(mouse2canine, HinC_res)
knitr::kable(rbind(CinC_enrichment_prob, HinC_enrichement_prob))

######################################
#
# (5) DUX4 targets gene set enrichment
#
######################################
HinH_res <- results(HinH.ens.dds, alpha=0.05, lfcThreshold=2)
HinH_up_reg <- as(HinH_res, "data.frame") %>%
  rownames_to_column(var="ENSEMBL") %>%
  dplyr::filter(padj < 0.05 & log2FoldChange >0)
seed_ensembl <- pull(HinH_up_reg, ENSEMBL)

dux4_homology <- human_inparanoid_homology(seed_ensembl)

CinC_enrichment_prob <- 
  geneset_hypergeometric(dux4_homology$human2canine, CinC_res)
HinC_enrichement_prob <- 
  geneset_hypergeometric(dux4_homology$human2canine, HinC_res)
knitr::kable(rbind(CinC_enrichment_prob, HinC_enrichement_prob))  



######################################
#
# (5) GO analysis
#
######################################
library(goseq)
library(ggplot2)
library(dplyr)
source(file.path(pkg_dir, "scripts", "tools.R"))

# (a) CinC
universe <- rownames(CinC.ens.dds)
selected <- CinC_res.2 %>% 
  dplyr::filter(CinC_padj < 0.05 & CinC_logFC > 0)  %>%
  pull(ENSEMBL) #1553
CinC_go <- .do_goseq(universe=universe, selected=selected, threshold_pval=0.01) %>%
  mutate(log10Pval = -10 * log10(over_represented_pvalue)) 
CinC_go_top35 <- CinC_go[1:35, ] %>%  
  dplyr::arrange(desc(over_represented_pvalue)) %>%
  mutate(term = factor(term, levels=term))

# where is the cutoff for fdr (adjusted p-value)?
fdr_n <- CinC_go_top35 %>% summarise(n = sum(fdr < 0.1))
gg <- ggplot(CinC_go_top35, aes(x=term, y=log10Pval)) +
  geom_bar(stat="identity") +
  geom_vline(xintercept = 35 - as.numeric(fdr_n) + 0.5, linetype="dashed", color="gray75") +
  coord_flip() +
  theme_minimal() +
  labs(y="-log10Pval", x="", title="CinC GO terms Enrichment")

pdf(file.path(fig_dir, "CinC_go.pdf"), width=12, height=5)
plot(gg)
dev.off()

# (b) HinC
universe <- rownames(HinC.ens.dds)
selected <- HinC_res.2 %>% 
  dplyr::filter(HinC_padj < 0.05 & HinC_logFC > 0)  %>%
  pull(ENSEMBL)
HinC_go <- .do_goseq(universe=universe, selected=selected, threshold_pval=0.01) %>%
  mutate(log10Pval = -10 * log10(over_represented_pvalue)) 
HinC_go_top35 <- HinC_go[1:35, ] %>% 
  dplyr::arrange(desc(over_represented_pvalue)) %>%
  mutate(term = factor(term, levels=term))

gg <- ggplot(HinC_go_top35, aes(x=term, y=log10Pval)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_minimal() +
  labs(y="-log10Pval", x="", title="HinC GO terms Enrichment")
pdf(file.path(fig_dir, "HinC_go.pdf"), width=12, height=5)
plot(gg)
dev.off()

# (c) HinC and CinC go combine
CinC_HinC_top35 <- CinC_go %>% top_n(35, log10Pval) %>%
  dplyr::arrange(log10Pval) %>%
  mutate(term = factor(term, levels=term)) %>%
  left_join(HinC_go, by="category", suffix=c("_CinC", "_HinC")) %>%
  dplyr::select(term_CinC, log10Pval_CinC, log10Pval_HinC) %>%
  rename(CinC=log10Pval_CinC, HinC=log10Pval_HinC) %>%
  gather(key="group", value="log10Pval", -term_CinC) 

gg <- ggplot(CinC_HinC_top35, aes(x=group, y=term_CinC)) +
  geom_point(aes(size=log10Pval)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  labs(title="Top GO in CinC")
pdf(file.path(fig_dir, "CinC_HinC_go_top35.pdf"))
plot(gg)
dev.off()  

HinC_CinC_top35 <- HinC_go %>% top_n(35, log10Pval) %>%
  dplyr::arrange(log10Pval) %>%
  mutate(term = factor(term, levels=term)) %>%
  left_join(CinC_go, by="category", suffix=c("_HinC", "_CinC")) %>%
  dplyr::select(term_HinC, log10Pval_CinC, log10Pval_HinC) %>%
  rename(CinC=log10Pval_CinC, HinC=log10Pval_HinC) %>%
  gather(key="group", value="log10Pval", -term_HinC) 

gg <- ggplot(HinC_CinC_top35, aes(x=group, y=term_HinC)) +
  geom_point(aes(size=log10Pval)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  labs(title="Top GO in HinC")
pdf(file.path(fig_dir, "HinC_CinC_go_top35.pdf"), width=4.8, height=7)
plot(gg)
dev.off()  

