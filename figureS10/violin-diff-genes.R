### getting packages

library(ggplot2)
library(ggsci)
library(data.table)
library(easyGgplot2)

### setting working directory

setwd("../data/")

### WT SMARCA4dtag

### reading data and merging

data_atac <- read.csv("differential_analysis.deseq_result.all_comparisons.csv")

data_want <- data_atac[data_atac$comparison_name %in% c("ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO"),]

data_want_sig <- data_want[data_want$padj < 0.01 & data_want$log2FoldChange < -1,]

data_chip <- read.table("macs2-log2-FC-rpm-based.txt",header = T)

data_chip_want <- data_chip[,c(1,28,31,34)]

data_atac_chip <- merge(data_want_sig,data_chip_want,"index")

data_gene <- read.csv("BAF_Timecourse_Set3.gene.annotation.singleline.txt")

colnames(data_gene)[1] <- "index"

data_rna <- read.table("rna-seq-all-comp-consensus-log2FConly.txt",header = T)

colnames(data_rna)[2] <- "gene_name"

data_rna_gene <- merge(data_gene,data_rna,"gene_name")

data_gene_diff <- read.table("WT_SMARCA4dtag_genes_differential.txt",header = F)

colnames(data_gene_diff) <- "gene_name"

data_gene_diff_unique <- data.frame(unique(data_gene_diff$gene_name))

colnames(data_gene_diff_unique) <- "gene_name"

data_merge_gene_rna_diff <- merge(data_gene_diff_unique,data_rna_gene,"gene_name")

data_rna_gene_want <- data_merge_gene_rna_diff[,c("gene_name","index","gene_id","log2FoldChange_WT_SMARCA4_F1_dTAG_47_2h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_3h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_6h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_24h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_72h_vs_HAP1_WT_SMARCA4_F1_DMSO")]

data_atac_chip_rna <- merge(data_atac_chip,data_rna_gene_want,"index")

### subsetting 

data_atac_chip_rna_down <- data_atac_chip_rna[data_atac_chip_rna$FC_WT_SMARCA4dTAG_DMSO_24h_H3K27ac < -1, ]

data_atac_chip_rna_down$condition <- "dtag-atac-decrease-24h-H3K27ac-chip-decrease"

print(length(data_atac_chip_rna_down$index))

data_atac_chip_rna_no_response <- data_atac_chip_rna[data_atac_chip_rna$FC_WT_SMARCA4dTAG_DMSO_24h_H3K27ac > -1 & data_atac_chip_rna$FC_WT_SMARCA4dTAG_DMSO_24h_H3K27ac < 1, ]

data_atac_chip_rna_no_response$condition <- "dtag-atac-decrease-24h-H3K27ac-chip-no-response"

print(length(data_atac_chip_rna_no_response$index))

data_plot_merge <- rbind(data_atac_chip_rna_down,data_atac_chip_rna_no_response)

### rna seq row wise

data_frame_list <- list(data_plot_merge[,c("log2FoldChange_WT_SMARCA4_F1_dTAG_47_2h_vs_HAP1_WT_SMARCA4_F1_DMSO","condition")],data_plot_merge[,c("log2FoldChange_WT_SMARCA4_F1_dTAG_47_3h_vs_HAP1_WT_SMARCA4_F1_DMSO","condition")],data_plot_merge[,c("log2FoldChange_WT_SMARCA4_F1_dTAG_47_6h_vs_HAP1_WT_SMARCA4_F1_DMSO","condition")],data_plot_merge[,c("log2FoldChange_WT_SMARCA4_F1_dTAG_47_24h_vs_HAP1_WT_SMARCA4_F1_DMSO","condition")],data_plot_merge[,c("log2FoldChange_WT_SMARCA4_F1_dTAG_47_72h_vs_HAP1_WT_SMARCA4_F1_DMSO","condition")])

data_plot <- rbindlist(data_frame_list,use.names=F)

colnames(data_plot)[1] <- "log2FC"

data_plot$time_point_rna_seq <- NA

data_plot$time_point_rna_seq[1:43] <- "2h"

data_plot$time_point_rna_seq[44:86] <- "3h"

data_plot$time_point_rna_seq[87:129] <- "6h"

data_plot$time_point_rna_seq[130:172] <- "24h"

data_plot$time_point_rna_seq[173:215] <- "72h"

data_plot_chose <- data_plot[data_plot$time_point_rna_seq %in% c("24h","72h"),]

### Making violin plots

level_vector <- c("24h","72h")

e <- ggplot(data_plot_chose, aes(x = factor(time_point_rna_seq,level_vector), y = log2FC)) +

 geom_violin(
  aes(color = condition), width = 0.3, size = 0.2,
  position = position_dodge(0.8)) + theme_bw() + ylim(-6.5,5) + stat_summary(fun=median, geom="point", size=2, color="black") +


geom_boxplot(
  aes(color = condition), width = 0.05, size = 0.2,
  position = position_dodge(0.8)) + theme_bw() + ylim(-6.5,5) + stat_summary(fun=median, geom="point", size=2, color="black") + theme(legend.position = "none")


ggsave("violin-WT-SMARCA4dtag-neg-log2FC1-atac-chip-diff-rna-seq-late-tp.pdf",height = 4, width = 5, dpi=300)
print(e)
dev.off()




### BRG BRMdtag

### reading data and merging

data_atac <- read.csv("differential_analysis.deseq_result.all_comparisons.csv")

data_want <- data_atac[data_atac$comparison_name %in% c("ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_72h_vs_DMSO"),]

data_want_sig <- data_want[data_want$padj < 0.01 & data_want$log2FoldChange < -1,]

data_chip <- read.table("/macs2-log2-FC-rpm-based.txt",header = T)

data_chip_want <- data_chip[,c(1,28,31,34)]

data_atac_chip <- merge(data_want_sig,data_chip_want,"index")

data_gene <- read.csv("BAF_Timecourse_Set3.gene.annotation.singleline.txt")

colnames(data_gene)[1] <- "index"

data_rna <- read.table("rna-seq-all-comp-consensus-log2FConly.txt",header = T)

colnames(data_rna)[2] <- "gene_name"

data_rna_gene <- merge(data_gene,data_rna,"gene_name")

data_gene_diff <- read.table("BRG_BRMdtag_genes_differential.txt",header = F)

colnames(data_gene_diff) <- "gene_name"

data_gene_diff_unique <- data.frame(unique(data_gene_diff$gene_name))

colnames(data_gene_diff_unique) <- "gene_name"

data_merge_gene_rna_diff <- merge(data_gene_diff_unique,data_rna_gene,"gene_name")

data_rna_gene_want <- data_merge_gene_rna_diff[,c("gene_name","index","gene_id","log2FoldChange_BRG_BRM_C8_dTAG_47_2h_vs_HAP1_BRG_BRM_C8_DMSO","log2FoldChange_BRG_BRM_C8_dTAG_47_3h_vs_HAP1_BRG_BRM_C8_DMSO","log2FoldChange_BRG_BRM_C8_dTAG_47_6h_vs_HAP1_BRG_BRM_C8_DMSO","log2FoldChange_BRG_BRM_C8_dTAG_47_24h_vs_HAP1_BRG_BRM_C8_DMSO","log2FoldChange_BRG_BRM_C8_dTAG_47_72h_vs_HAP1_BRG_BRM_C8_DMSO")]

data_atac_chip_rna <- merge(data_atac_chip,data_rna_gene_want,"index")

### subsetting 

data_atac_chip_rna_down <- data_atac_chip_rna[data_atac_chip_rna$FC_SMARCA4KO_SMARCA2dTAG_DMSO_24h_H3K27ac < -1, ]

data_atac_chip_rna_down$condition <- "dtag-atac-decrease-24h-H3K27ac-chip-decrease"

print(length(data_atac_chip_rna_down$index))

data_atac_chip_rna_no_response <- data_atac_chip_rna[data_atac_chip_rna$FC_SMARCA4KO_SMARCA2dTAG_DMSO_24h_H3K27ac > -1 & data_atac_chip_rna$FC_SMARCA4KO_SMARCA2dTAG_DMSO_24h_H3K27ac < 1, ]

data_atac_chip_rna_no_response$condition <- "dtag-atac-decrease-24h-H3K27ac-chip-no-response"

print(length(data_atac_chip_rna_no_response$index))

data_plot_merge <- rbind(data_atac_chip_rna_down,data_atac_chip_rna_no_response)

### rna seq row wise

data_frame_list <- list(data_plot_merge[,c("log2FoldChange_BRG_BRM_C8_dTAG_47_2h_vs_HAP1_BRG_BRM_C8_DMSO","condition")],data_plot_merge[,c("log2FoldChange_BRG_BRM_C8_dTAG_47_3h_vs_HAP1_BRG_BRM_C8_DMSO","condition")],data_plot_merge[,c("log2FoldChange_BRG_BRM_C8_dTAG_47_6h_vs_HAP1_BRG_BRM_C8_DMSO","condition")],data_plot_merge[,c("log2FoldChange_BRG_BRM_C8_dTAG_47_24h_vs_HAP1_BRG_BRM_C8_DMSO","condition")],data_plot_merge[,c("log2FoldChange_BRG_BRM_C8_dTAG_47_72h_vs_HAP1_BRG_BRM_C8_DMSO","condition")])

data_plot <- rbindlist(data_frame_list,use.names=F)

colnames(data_plot)[1] <- "log2FC"

data_plot$time_point_rna_seq <- NA

data_plot$time_point_rna_seq[1:43] <- "2h"

data_plot$time_point_rna_seq[44:86] <- "3h"

data_plot$time_point_rna_seq[87:129] <- "6h"

data_plot$time_point_rna_seq[130:172] <- "24h"

data_plot$time_point_rna_seq[173:215] <- "72h"

data_plot_chose <- data_plot[data_plot$time_point_rna_seq %in% c("24h","72h"),]

level_vector <- c("24h","72h")

### Making violin plots

e <- ggplot(data_plot_chose, aes(x = factor(time_point_rna_seq,level_vector), y = log2FC)) +
  
  geom_violin(
    aes(color = condition), width = 0.3, size = 0.2,
    position = position_dodge(0.8)) + theme_bw() + ylim(-6.5,5) +
  
  geom_boxplot(
    aes(color = condition), width = 0.05, size = 0.2,
    position = position_dodge(0.8)) + theme_bw() + ylim(-6.5,5) + theme(legend.position = "none")

ggsave("violin-BRG-BRMdtag-neg-log2FC1-atac-chip-diff-rna-seq-late-tp.pdf",height = 4, width = 5, dpi=300)
print(e)
dev.off()








### CC1 CC2dtag

### reading data and merging

data_atac <- read.csv("differential_analysis.deseq_result.all_comparisons.csv")

data_want <- data_atac[data_atac$comparison_name %in% c("ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_3h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_6h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_24h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_72h_vs_combined_control"),]

data_want_sig <- data_want[data_want$padj < 0.01 & data_want$log2FoldChange < -1,]

data_chip <- read.table("macs2-log2-FC-rpm-based.txt",header = T)

data_chip_want <- data_chip[,c(1,28,31,34)]

data_atac_chip <- merge(data_want_sig,data_chip_want,"index")

data_gene <- read.csv("BAF_Timecourse_Set3.gene.annotation.singleline.txt")

colnames(data_gene)[1] <- "index"

data_rna <- read.table("rna-seq-all-comp-consensus-log2FConly.txt",header = T)

colnames(data_rna)[2] <- "gene_name"

data_rna_gene <- merge(data_gene,data_rna,"gene_name")

data_gene_diff <- read.table("CC1_CC2dtag_genes_differential.txt",header = F)

colnames(data_gene_diff) <- "gene_name"

data_gene_diff_unique <- data.frame(unique(data_gene_diff$gene_name))

colnames(data_gene_diff_unique) <- "gene_name"

data_merge_gene_rna_diff <- merge(data_gene_diff_unique,data_rna_gene,"gene_name")

data_rna_gene_want <- data_merge_gene_rna_diff[,c("gene_name","index","gene_id","log2FoldChange_CC1_CC2_F1_dTAG_47_2h_vs_HAP1_CC1_CC2_F1_DMSO","log2FoldChange_CC1_CC2_F1_dTAG_47_3h_vs_HAP1_CC1_CC2_F1_DMSO","log2FoldChange_CC1_CC2_F1_dTAG_47_6h_vs_HAP1_CC1_CC2_F1_DMSO","log2FoldChange_CC1_CC2_F1_dTAG_47_24h_vs_HAP1_CC1_CC2_F1_DMSO","log2FoldChange_CC1_CC2_F1_dTAG_47_72h_vs_HAP1_CC1_CC2_F1_DMSO")]

data_atac_chip_rna <- merge(data_atac_chip,data_rna_gene_want,"index")

### subsetting 

data_atac_chip_rna_down <- data_atac_chip_rna[data_atac_chip_rna$FC_SMARCC1KO_SMARCC2dTAG_DMSO_24h_H3K27ac < -1, ]

data_atac_chip_rna_down$condition <- "dtag-atac-decrease-24h-H3K27ac-chip-decrease"

print(length(data_atac_chip_rna_down$index))

data_atac_chip_rna_no_response <- data_atac_chip_rna[data_atac_chip_rna$FC_SMARCC1KO_SMARCC2dTAG_DMSO_24h_H3K27ac > -1 & data_atac_chip_rna$FC_SMARCC1KO_SMARCC2dTAG_DMSO_24h_H3K27ac < 1, ]

data_atac_chip_rna_no_response$condition <- "dtag-atac-decrease-24h-H3K27ac-chip-no-response"

print(length(data_atac_chip_rna_no_response$index))

data_plot_merge <- rbind(data_atac_chip_rna_down,data_atac_chip_rna_no_response)

### rna seq row wise

data_frame_list <- list(data_plot_merge[,c("log2FoldChange_CC1_CC2_F1_dTAG_47_2h_vs_HAP1_CC1_CC2_F1_DMSO","condition")],data_plot_merge[,c("log2FoldChange_CC1_CC2_F1_dTAG_47_3h_vs_HAP1_CC1_CC2_F1_DMSO","condition")],data_plot_merge[,c("log2FoldChange_CC1_CC2_F1_dTAG_47_6h_vs_HAP1_CC1_CC2_F1_DMSO","condition")],data_plot_merge[,c("log2FoldChange_CC1_CC2_F1_dTAG_47_24h_vs_HAP1_CC1_CC2_F1_DMSO","condition")],data_plot_merge[,c("log2FoldChange_CC1_CC2_F1_dTAG_47_72h_vs_HAP1_CC1_CC2_F1_DMSO","condition")])

data_plot <- rbindlist(data_frame_list,use.names=F)

colnames(data_plot)[1] <- "log2FC"

data_plot$time_point_rna_seq <- NA

data_plot$time_point_rna_seq[1:203] <- "2h"

data_plot$time_point_rna_seq[204:406] <- "3h"

data_plot$time_point_rna_seq[407:609] <- "6h"

data_plot$time_point_rna_seq[610:812] <- "24h"

data_plot$time_point_rna_seq[813:1015] <- "72h"

data_plot_chose <- data_plot[data_plot$time_point_rna_seq %in% c("24h","72h"),]

level_vector <- c("24h","72h")

### Making violin plots


e <- ggplot(data_plot_chose, aes(x = factor(time_point_rna_seq,level_vector), y = log2FC)) +
  
  geom_violin(
    aes(color = condition), width = 0.3, size = 0.2,
    position = position_dodge(0.8)) + theme_bw() + ylim(-6.5,5) +
  
  geom_boxplot(
    aes(color = condition), width = 0.05, size = 0.2,
    position = position_dodge(0.8)) + theme_bw() + ylim(-6.5,5) + theme(legend.position = "none")

ggsave("violin-CC1-CC2dtag-neg-log2FC1-atac-chip-diff-rna-seq-late-tp.pdf",height = 4, width = 5, dpi=300)
print(e)
dev.off()