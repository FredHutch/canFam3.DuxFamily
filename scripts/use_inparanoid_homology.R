#' this script can later be rename: 05_inparanoid_homolgy
#' This script use Bioconductor's homology package built by inparanoid algorithm
#' We need human to canine


library(DBI)
library("hom.Hs.inp.db")
library(dplyr)
library(tidyr)
library(purrr)

library(lattice)
library(latticeExtra)

library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Cf.eg.db)

# "hom.Hs.inpCANFA: human to canine
pkg_dir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/canFam3.DuxFamily"
#local_dir <- "/Users/cwon2/CompBio/canFam3.DuxFamily"
data_dir <- file.path(pkg_dir, "data")
source(file.path(pkg_dir, "scripts", "inparanoid_homology.R"))
source(file.path(pkg_dir, "scripts", "tools.R"))

# loading dataset
load(file.path(data_dir,"C2C12.ens.ddsl2.rda"))
load(file.path(data_dir, "CinC.ens.dds.rda"))
rownames(CinC.ens.dds) <- sapply(strsplit(rownames(CinC.ens.dds), fixed=TRUE, "."), "[", 1)
load(file.path(data_dir, "HinC.ens.dds.rda"))
rownames(HinC.ens.dds) <- sapply(strsplit(rownames(HinC.ens.dds), fixed=TRUE, "."), "[", 1)
load(file.path(data_dir, "HinH.ens.dds.rda"))
load(file.path(data_dir, "human_cleavage.rda"))
my_color <- c('#999999','#E69F00', '#56B4E9')

#
# cleavage_homology and 2C-lie
#
cleavage_homology <- human_inparanoid_homology(as.character(human_cleavage$GeneID))
# how many of them are differentially expressed? run GSEA
# like2c_homology

#
# whole gene set homology
#
load(file.path(data_dir, "HinH.ens.SE.rda"))
no_keep <- grepl("LRG", rownames(HinH.ens.SE))
human_universe <- rownames(HinH.ens.SE)[!no_keep]
universe_homology <- human_inparanoid_homology(human_universe)


#
# get DESeq2 resuls for CinC, HinC, HinH, MinM
#
library(DESeq2)
library(tidyverse)
HinH_res <- as(results(HinH.ens.dds, lfcThreshold=1, alpha=0.05), "data.frame") %>%
  rownames_to_column(var="HUMAN_ENSEMBL") %>%
  dplyr::select(HUMAN_ENSEMBL, log2FoldChange, padj) %>%
  rename(HinH_logFC=log2FoldChange, HinH_padj=padj) %>%
  dplyr::mutate(HinH_DE=HinH_padj < 0.05)

res <- results(CinC.ens.dds, lfcThreshold=1, alpha=0.05)
summary(res)
CinC_res <- as(res, "data.frame") %>%
  rownames_to_column(var="CANINE_ENSEMBL") %>%
  dplyr::select(CANINE_ENSEMBL, log2FoldChange, padj) %>%
  rename(CinC_logFC=log2FoldChange, CinC_padj=padj) %>%
  dplyr::mutate(CinC_DE=CinC_padj < 0.05)

res <- results(HinC.ens.dds, lfcThreshold=1, alpha=0.05)
summary(res)
HinC_res <- as(res, "data.frame") %>%
  rownames_to_column(var="CANINE_ENSEMBL") %>%
  dplyr::select(CANINE_ENSEMBL, log2FoldChange, padj) %>%
  rename(HinC_logFC=log2FoldChange, HinC_padj=padj) %>%
  dplyr::mutate(HinC_DE=HinC_padj < 0.05)


MinM_res <- as(results(C2C12.ens.ddsl2[[2]], lfcThreshold=1, alpha=0.05), "data.frame") %>%
  rownames_to_column(var="MOUSE_ENSEMBL") %>%
  dplyr::select(MOUSE_ENSEMBL, log2FoldChange, padj) %>%
  rename(MinM_logFC=log2FoldChange, MinM_padj=padj) %>%
  dplyr::mutate(MinM_DE=MinM_padj < 0.05)

#  sanity check
sum(HinH_res$HUMAN_ENSEMBL %in% universe_homology[[1]]$HUMAN_ENSEMBL)
sum(CinC_res$CANINE_ENSEMBL %in% universe_homology[[1]]$CANINE_ENSEMBL)

#
# cross-species homology
#


#
# (1) HinC and CinC on human homology set
#
HinC_CinC <- HinC_res %>%
  left_join(CinC_res, by="CANINE_ENSEMBL") %>%
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

#
#  (2) HinH and CinC on human homology set
#
HinH_CinC <- universe_homology[["human2canine"]] %>%
  inner_join(HinH_res, by="HUMAN_ENSEMBL") %>% #1685
  left_join(CinC_res, by="CANINE_ENSEMBL") %>%# 1685
  dplyr::mutate(both_DE= HinH_DE & CinC_DE) %>%
  dplyr::mutate(DE_status = "neither") %>%
  dplyr::mutate(DE_status = replace(DE_status, HinH_DE, "either")) %>%
  dplyr::mutate(DE_status = replace(DE_status, CinC_DE, "either")) %>%
  dplyr::mutate(DE_status = replace(DE_status, both_DE, "both")) %>%
  dplyr::mutate(DE_status = factor(DE_status, levels=c("neither", "either", "both")))
table(HinH_CinC$DE_status)

gg <- ggplot(HinH_CinC, aes(x=CinC_logFC, y=HinH_logFC)) +
  geom_point(size=1, alpha=0.4, aes(color=DE_status)) +
  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
  geom_smooth(method="lm", se=FALSE, show.legend=FALSE, color="gray70", alpha=0.5,
              linetype="dashed") +
  labs(title="HinH vs CinC logFC", color="DE status") +
  theme_bw()+
  theme(legend.position = c(0.15, 0.8)) 
df <- HinH_CinC %>% dplyr::select(HinH_logFC, CinC_logFC) %>%
  rename(y=HinH_logFC, x=CinC_logFC)  
gg <- gg + geom_text(x = 1, y = 7.5, vjust=0, hjust=0, 
  label = lm_eqn(df), parse = TRUE, color="gray50") 

pdf(file.path(pkg_dir, "figures", "HinH_CinC_logFC_scatter.pdf"),
    width=4, height=4)  
plot(gg)
dev.off()

#
# (2) HinH and MinM on human homology set
#
HinH_MinM <- universe_homology[[2]] %>%
  inner_join(HinH_res, by="HUMAN_ENSEMBL") %>% #1208
  left_join(MinM_res, by="MOUSE_ENSEMBL") %>%
  dplyr::mutate(both_DE= HinH_DE & MinM_DE) %>%
  dplyr::mutate(DE_status = "neither") %>%
  dplyr::mutate(DE_status = replace(DE_status, HinH_DE, "either")) %>%
  dplyr::mutate(DE_status = replace(DE_status, MinM_DE, "either")) %>%
  dplyr::mutate(DE_status = replace(DE_status, both_DE, "both")) %>%
  dplyr::mutate(DE_status = factor(DE_status, levels=c("neither", "either", "both")))
table(HinH_MinM$DE_status)

gg <- ggplot(HinH_MinM, aes(x=MinM_logFC, y=HinH_logFC)) +
  geom_point(size=1, alpha=0.4, aes(color=DE_status)) +
  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
  geom_smooth(method="lm", se=FALSE, show.legend=FALSE, color="gray70", alpha=0.5,
              linetype="dashed") +
  labs(title="HinH vs MinM logFC", color="DE status") +
  theme_bw()+
  theme(legend.position = c(0.15, 0.8)) 
df <- HinH_MinM %>% dplyr::select(HinH_logFC, MinM_logFC) %>%
  rename(y=HinH_logFC, x=MinM_logFC)  
gg <- gg + geom_text(x = 1, y = 7.5, vjust=0, hjust=0, 
  label = lm_eqn(df), parse = TRUE, color="gray50")   

pdf(file.path(pkg_dir, "figures", "HinH_MinM_logFC_scatter.pdf"),
    width=4, height=4)  
plot(gg)
dev.off()

#
# (3) HinH and HinC on human homology set
#
HinH_HinC <- universe_homology[["human2canine"]] %>%
  inner_join(HinH_res, by="HUMAN_ENSEMBL") %>% 
  left_join(HinC_res, by="CANINE_ENSEMBL") %>%
  dplyr::mutate(both_DE= HinH_DE & HinC_DE) %>%
  dplyr::mutate(DE_status = "neither") %>%
  dplyr::mutate(DE_status = replace(DE_status, HinH_DE, "either")) %>%
  dplyr::mutate(DE_status = replace(DE_status, HinC_DE, "either")) %>%
  dplyr::mutate(DE_status = replace(DE_status, both_DE, "both")) %>%
  dplyr::mutate(DE_status = factor(DE_status, levels=c("neither", "either", "both")))
table(HinH_HinC$DE_status)

gg <- ggplot(HinH_HinC, aes(x=HinC_logFC, y=HinH_logFC)) +
  geom_point(size=1, alpha=0.4, aes(color=DE_status)) +
  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
  geom_smooth(method="lm", se=FALSE, show.legend=FALSE, color="gray70", alpha=0.5,
              linetype="dashed") +
  labs(title="HinH vs HinC logFC", color="DE status") +
  theme_bw()+
  theme(legend.position = c(0.15, 0.8)) 
df <- HinH_HinC %>% dplyr::select(HinH_logFC, HinC_logFC) %>%
  rename(y=HinH_logFC, x=HinC_logFC)  
gg <- gg + geom_text(x = 1, y = 7.5, vjust=0, hjust=0, 
  label = lm_eqn(df), parse = TRUE, color="gray50")   

pdf(file.path(pkg_dir, "figures", "HinH_HinC_logFC_scatter.pdf"),
    width=4, height=4)  
plot(gg)
dev.off()


#
# cleavage-stage gene signature only
#



