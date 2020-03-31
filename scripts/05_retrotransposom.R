# This script performs repeat element analysis

pkg_dir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/canFam3.DuxFamily"
rmsk_dir <- "/fh/fast/tapscott_s/CompBio/R_packages/rmskStats"
fig_dir <- file.path(pkg_dir, "figures")
data_dir <- file.path(pkg_dir, "data")

library(DESeq2)
library(tidyverse)
library(tidyr)
library(ggrepel)

load(file.path(pkg_dir, "data", "rmsk.dds.rda"))
load(file.path(pkg_dir, "data", "rmsk.res.rda"))

# 
# (1)DESeq2 results
#


# (1a) MA plot
lapply(names(rmsk.dds[1:2]), function(name) {
    res <- results(rmsk.dds[[name]], alpha=0.05, lfcThreshold = 0.5)
    summary <- as(res, "data.frame") %>%
      summarise(up = sum(padj < 0.05 & log2FoldChange > 0.95, na.rm=TRUE), 
      down = sum(padj < 0.05 & log2FoldChange < 0.95, na.rm=TRUE))

    file_name <- paste0("rmsk_", name, "_MAplot.pdf")
    y_text <- c(down=-4, up=4)

    if (name == "CinC") y_text <- c(down=-2, up=2)

    pdf(file.path(fig_dir, file_name), height=4, width=4.5)
    plotMA(res, alpha = 0.05, 
           main=paste0(name, ": ", nrow(res), " annotated elements"))
    label <- sprintf("LFC > 0.5 (up) : %d, %1.1f%%", 
                     pull(summary, up), pull(summary, up)/nrow(res) *100)           
    text(x=80, y=y_text["up"],   labels=label, adj=c(0, 0))
    label <- sprintf("LFC < 0.5 (down) : %d, %1.1f%%", 
                     pull(summary, down), pull(summary, down)/nrow(res) *100)     
    text(x=80, y=y_text["down"], labels=label, adj=c(0, 1))
    dev.off()  
})

# (1b) Enrichment Analysis

source(file.path(pkg_dir, "scripts", "tools.R")) # .rmsk_enrichment
HinC_rmsk <- .rmsk_enrichment(rmsk.dds[["HinC"]])
CinC_rmsk <- .rmsk_enrichment(rmsk.dds[["CinC"]])

# I would like to know the enrichment analysis results differ by different
# threshold? not much different (.rmsk_enrichement_alt()) Stable 
# with differnt thresholds

# satellite: 
HinC_rmsk$res_df %>% dplyr::filter(repFamily == "Satellite")
CinC_rmsk$res_df %>% dplyr::filter(repFamily == "Satellite")

# (1c) scatter plot of HinC and CinC; highlight LTR: ERVL-MaLR
comb <- inner_join(HinC_rmsk$res_df, CinC_rmsk$res_df, 
                   by="repName", suffix=c("_HinC", "_CinC")) %>%
  dplyr::mutate(show_name = as.character(repFamily_HinC)) %>%
  dplyr::mutate(show_name = ifelse(abs(log2FoldChange_HinC-log2FoldChange_CinC) < 1.5, 
                                   "", show_name)) %>%      
  dplyr::mutate(DE_status = "neither") %>%
  dplyr::mutate(DE_status = replace(DE_status, padj_HinC < 0.05 | padj_CinC < 0.05, "either")) %>%
  dplyr::mutate(DE_status = replace(DE_status, padj_HinC < 0.05 & padj_CinC < 0.05, "both")) %>%
  dplyr::mutate(DE_status = factor(DE_status, levels=c("neither", "either", "both")))                                          

comb%>% summarise(both=sum(DE_status==both, na.rm=TRUE)/n())

comb %>% dplyr::filter(log2FoldChange_HinC - log2FoldChange_CinC > 1.5)  
# label the |lfc| > 1.5                
library(ggrepel)

gg1 <- ggplot(comb, aes(x=log2FoldChange_CinC, y=log2FoldChange_HinC)) +
  geom_point(size=1, alpha=0.7, aes(color=DE_status)) +
  scale_color_manual(values=c("#999999","#E69F00", "#56B4E9")) +
  geom_smooth(method="lm", se=FALSE, show.legend=FALSE, color="gray66", alpha=0.5,
              linetype="dashed") +
  geom_abline(intercept=0, slope=1, color="gray75", alpha=0.2, size=0.5) +              
  labs(title="RMSK: HinC vs CinC log2FC") +
  theme_bw()+
  theme(legend.position = c(0.1, 0.88)) +
  geom_text_repel(aes(label=show_name), size=2.5, show.legend=FALSE)

df <- comb %>% dplyr::select(log2FoldChange_HinC, log2FoldChange_CinC) %>%
  rename(y=log2FoldChange_HinC, x=log2FoldChange_CinC)  
gg1 <- gg1 + geom_text(x = 7, y = 8.3, vjust=0, hjust=0.5, 
  label = lm_eqn(df), parse = TRUE, color="gray65", alpha=0.5)

pdf(file.path(fig_dir, "rmsk_HinC_CinC_log2FC_scatter.pdf"),
    width=6, height=6)  
gg1
dev.off()

#
# (2) enrichment comparison
#

# (2a) over-represented repeat families and classes
family_res <- full_join(HinC_rmsk$repFamily_summary, CinC_rmsk$repFamily_summary, 
                        by="repFamily", suffix=c("_HinC", "_CinC")) %>%
  dplyr::select(repFamily, enrichment_prob_HinC, enrichment_prob_CinC,
    up_HinC, up_CinC, enrichment_mu_HinC, enrichment_mu_CinC) %>% 
  tidyr::gather(key="sample", value="pval", -repFamily, 
                -up_HinC, -up_CinC, -enrichment_mu_HinC, -enrichment_mu_CinC) %>%
  dplyr::mutate(sample = ifelse(grepl("HinC", sample), "HinC", "CinC")) %>%
  dplyr::mutate(log10Pval = -10*log10(pval)) %>%
  dplyr::mutate(up = ifelse(sample=="HinC", up_HinC, up_CinC),
                mu = ifelse(sample=="HinC", enrichment_mu_HinC, enrichment_mu_CinC)) %>%
  dplyr::mutate(overrepresented = pval < 0.1 & up > mu) 

gg <- ggplot(family_res, aes(x=repFamily, y=sample)) +
  geom_point(aes(size=log10Pval, color=overrepresented, fill=overrepresented)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="bottom") +
  labs(title=expression("Repeat Family GSEA: " ~ -10 * Log[10] ~ "(p-value)"),
       size = expression(-10 * Log[10] ~ "(p-value)")) +
  scale_fill_manual(name ="over-represented",  values=c("ivory4", "firebrick4")) +
  scale_color_manual(name="over-represented", values=c("ivory4", "firebrick4"))

pdf(file.path(fig_dir, "rmsk_repFamily_log10Pval.2.pdf"), width=8, height=3.8)
gg
dev.off()

class_res <- full_join(HinC_rmsk$repClass_summary, CinC_rmsk$repClass_summary, 
                        by="repClass", suffix=c("_HinC", "_CinC")) %>%
  dplyr::select(repClass, enrichment_prob_HinC, enrichment_prob_CinC,
    up_HinC, up_CinC, enrichment_mu_HinC, enrichment_mu_CinC) %>% 
  tidyr::gather(key="sample", value="pval", -repClass, 
                -up_HinC, -up_CinC, -enrichment_mu_HinC, -enrichment_mu_CinC) %>%
  dplyr::mutate(sample = ifelse(grepl("HinC", sample), "HinC", "CinC")) %>%
  dplyr::mutate(log10Pval = -10*log10(pval)) %>%
  dplyr::mutate(up = ifelse(sample=="HinC", up_HinC, up_CinC),
                mu = ifelse(sample=="HinC", enrichment_mu_HinC, enrichment_mu_CinC)) %>%
  dplyr::mutate(overrepresented = pval < 0.1 & up > mu)                         

gg <- ggplot(class_res, aes(x=repClass, y=sample)) +
  geom_point(aes(size=log10Pval, color=overrepresented, fill=overrepresented), show.legend=FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="bottom") +
  labs(title=expression("Repeat Class GSEA: " ~ -10 * Log[10] ~ "(p-value)"),
       size = expression(-10 * Log[10] ~ "(p-value)")) +
  scale_fill_manual(name ="over-represented",  values=c("ivory4", "firebrick4")) +
  scale_color_manual(name="over-represented", values=c("ivory4", "firebrick4"))

pdf(file.path(fig_dir, "rmsk_repClass_log10Pval.pdf"), width=4.5, height=2)
gg
dev.off()

# (2b) visualization of LTR families in terms of log2FC
LTR <- comb %>% dplyr::filter(repClass_HinC == "LTR") %>%
  dplyr::mutate(repFamily_HinC = factor(repFamily_HinC)) %>%
  dplyr::select(repName, repFamily_HinC, repClass_HinC, log2FoldChange_HinC, log2FoldChange_CinC, 
                padj_HinC, padj_CinC) %>%
  dplyr::rename(repFamily=repFamily_HinC, repClass=repClass_HinC) %>%
  tidyr::gather(key="sample", value="log2FC", -repName, -repFamily, -repClass, -padj_HinC, -padj_CinC) %>%
  dplyr::mutate(sample = ifelse(grepl("HinC", sample), "HinC", "CinC")) %>%
  dplyr::mutate(padj = ifelse(sample == "HinC", padj_HinC, padj_CinC)) %>%
  dplyr::mutate(DE_status = padj < 0.05) %>% # fix repFamily: put *? to others
  dplyr::mutate(repFamily=as.character(repFamily)) %>%
  dplyr::mutate(repFamily = 
                replace(repFamily, which(repFamily %in% c("ERV1?", "ERVL?", "Gypsy?", "LTR")), "Others"))

gg <- ggplot(LTR, aes(x=repFamily, y=log2FC)) +
  #geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_violin(width=0.7) +
  #geom_jitter(size=0.6, alpha=0.7, aes(color=DE_status), show.legend=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, 
               aes(color=DE_status, fill=DE_status), 
               alpha=0.9, show.legend=FALSE) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=2, color="red") +               
  theme_bw() +
  facet_wrap( ~ sample, nrow=2) +
  scale_fill_manual(values=c("ivory4", "firebrick4")) +
  scale_color_manual(values=c("ivory4", "firebrick4"))  +
  labs(title="Repeat families in LTR class", x="family")

pdf(file.path(fig_dir, "rmsk_repFamily_boxplot.pdf"), width=4, height=6)
gg
dev.off()

#
# binding sites in relation to LTR family
#
chip_dir <- "/fh/fast/tapscott_s/CompBio/ChIP-Seq/canFam3.ImmCnMb.dux"
load(file.path(chip_dir, "data", "peaks_list.rda"))

# ? how many peaks/percentage are overlapping with LTR 
sapply(peaks_list, function(rng) {
    n <- sum(rng$repClass == "LTR", na.rm=TRUE) 
    return(c(total=length(rng), n_LTR=n, freq_LTR=n/length(rng)))
})

# tidy up frequency: CinC stronger in ERVL; HinC stronger in ERVL-MaLR
LTR_freq <- sapply(peaks_list[-3], function(rng) {
    100 * table(rng$repFamily)[c("ERV1", "ERVL", "ERVL-MaLR")] / length(rng)
})
LTR_freq <- as.data.frame(t(LTR_freq)) %>%
  rownames_to_column(var="sample_name") %>%
  dplyr::mutate(trans_factor = sapply(strsplit(sample_name, "_", fixed=TRUE), "[[", 3)) %>%
  dplyr::mutate(clone_type = sapply(strsplit(sample_name, "_", fixed=TRUE), "[[", 4)) %>%
  tidyr::gather(key="Family", value="Frequency", -trans_factor, -clone_type, -sample_name) %>%
  dplyr::mutate(Sample = paste0(trans_factor, "_", clone_type)) %>%
  dplyr::mutate(Sample=factor(Sample, 
                            levels = c("hDUX4_poly", "hDUX4_mono", "CCH_poly", "CCH_mono", "CALTh_poly")))
gg <- ggplot(LTR_freq, aes(x=Sample, y=Frequency, fill=trans_factor)) +
  geom_bar(stat="identity", width=0.8) +
  facet_wrap( ~ Family, nrow=1) +
  geom_text(aes(label=format(Frequency, digits=2)), vjust=1.6, color="gray25", size=3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  labs(y="Frequency (%)") +       
  scale_fill_manual(values=c("#0072B2", "#009E73")) 
pdf(file.path(fig_dir, "rmsk_LTR_peaks_freq_bar.pdf"), width=6, height=3)
gg
dev.off()  

# tidy up average fold enrichment

#
# other threshold
#
as(results(rmsk.dds[["HinC"]], alpha=0.05), "data.frame") %>%
  summarise(up = sum(padj < 0.05 & log2FoldChange > 0.95, na.rm=TRUE), 
            down = sum(padj < 0.05 & log2FoldChange < 0.95, na.rm=TRUE))

as(results(rmsk.dds[["CinC"]], alpha=0.05), "data.frame") %>%
  summarise(up = sum(padj < 0.05 & log2FoldChange > 0.95, na.rm=TRUE), 
            down = sum(padj < 0.05 & log2FoldChange < 0.95, na.rm=TRUE))

