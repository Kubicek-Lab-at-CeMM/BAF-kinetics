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

data_sup_enh <- read.table("superenhancer_consensus_overlap.bed",header = F,sep = "\t")

colnames(data_sup_enh) <- c("chr","pos1","pos2")

data_sup_enh$index <- paste(data_sup_enh$chr,":",data_sup_enh$pos1,"-",data_sup_enh$pos2,sep = "")

data_atac_sup <- merge(data_atac,data_sup_enh,"index")

data_want <- data_atac_sup[data_atac_sup$comparison_name %in% c("ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO"),]

data_want_sig <- data_want[data_want$padj < 0.01 & data_want$log2FoldChange < -1,]

data_chip <- read.table("macs2-log2-FC-rpm-based.txt",header = T)

data_chip_want <- data_chip[,c(1,28,31,34)]

data_atac_chip <- merge(data_want_sig,data_chip_want,"index")

data_gene <- read.csv("BAF_Timecourse_Set3.gene.annotation.singleline.txt")

colnames(data_gene)[1] <- "index"

data_rna <- read.table("rna-seq-all-comp-consensus-log2FConly.txt",header = T)

colnames(data_rna)[2] <- "gene_name"

data_rna_gene <- merge(data_gene,data_rna,"gene_name")

data_rna_gene_want <- data_rna_gene[,c("gene_name","index","gene_id","log2FoldChange_WT_SMARCA4_F1_dTAG_47_2h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_3h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_6h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_24h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_72h_vs_HAP1_WT_SMARCA4_F1_DMSO")]

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

data_plot$time_point_rna_seq[1:37] <- "2h"

data_plot$time_point_rna_seq[38:74] <- "3h"

data_plot$time_point_rna_seq[75:111] <- "6h"

data_plot$time_point_rna_seq[112:148] <- "24h"

data_plot$time_point_rna_seq[149:185] <- "72h"

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"]))))

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"]))))

wilcox.test(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"])),unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"])))

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"]))))

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"]))))

wilcox.test(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"])),unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"])))


### Making violin plots

level_vector <- c("24h","72h")

e <- ggplot(data_plot[data_plot$time_point_rna_seq %in% level_vector,], aes(x = factor(time_point_rna_seq,level_vector), y = log2FC)) +

 geom_violin(
  aes(color = condition), width = 0.3, size = 0.2,
  position = position_dodge(0.8)) + theme_bw() + ylim(-3,1) +
  
  geom_boxplot(
    aes(color = condition), width = 0.05, size = 0.2,
    position = position_dodge(0.8)) + theme_bw() + ylim(-3,1) + theme(legend.position = "none")

ggsave("violin-WT-SMARCA4dtag-neg-log2FC1-atac-chip-all-rna-seq-superenhancer.pdf",height = 4, width = 5, dpi=300)
print(e)
dev.off()




### BRG BRMdtag

### reading data and merging

data_atac <- read.csv("differential_analysis.deseq_result.all_comparisons.csv")

data_sup_enh <- read.table("superenhancer_consensus_overlap.bed",header = F,sep = "\t")

colnames(data_sup_enh) <- c("chr","pos1","pos2")

data_sup_enh$index <- paste(data_sup_enh$chr,":",data_sup_enh$pos1,"-",data_sup_enh$pos2,sep = "")

data_atac_sup <- merge(data_atac,data_sup_enh,"index")

data_want <- data_atac_sup[data_atac_sup$comparison_name %in% c("ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_72h_vs_DMSO"),]

data_want_sig <- data_want[data_want$padj < 0.01 & data_want$log2FoldChange < -1,]

data_chip <- read.table("macs2-log2-FC-rpm-based.txt",header = T)

data_chip_want <- data_chip[,c(1,28,31,34)]

data_atac_chip <- merge(data_want_sig,data_chip_want,"index")

data_gene <- read.csv("BAF_Timecourse_Set3.gene.annotation.singleline.txt")

colnames(data_gene)[1] <- "index"

data_rna <- read.table("rna-seq-all-comp-consensus-log2FConly.txt",header = T)

colnames(data_rna)[2] <- "gene_name"

data_rna_gene <- merge(data_gene,data_rna,"gene_name")

data_rna_gene_want <- data_rna_gene[,c("gene_name","index","gene_id","log2FoldChange_BRG_BRM_C8_dTAG_47_2h_vs_HAP1_BRG_BRM_C8_DMSO","log2FoldChange_BRG_BRM_C8_dTAG_47_3h_vs_HAP1_BRG_BRM_C8_DMSO","log2FoldChange_BRG_BRM_C8_dTAG_47_6h_vs_HAP1_BRG_BRM_C8_DMSO","log2FoldChange_BRG_BRM_C8_dTAG_47_24h_vs_HAP1_BRG_BRM_C8_DMSO","log2FoldChange_BRG_BRM_C8_dTAG_47_72h_vs_HAP1_BRG_BRM_C8_DMSO")]

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

data_plot$time_point_rna_seq[1:101] <- "2h"

data_plot$time_point_rna_seq[102:202] <- "3h"

data_plot$time_point_rna_seq[203:303] <- "6h"

data_plot$time_point_rna_seq[304:404] <- "24h"

data_plot$time_point_rna_seq[405:505] <- "72h"

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"]))))

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"]))))

wilcox.test(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"])),unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"])))

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"]))))

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"]))))

wilcox.test(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"])),unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"])))


### Making violin plots

#level_vector <- c("2h","3h","6h","24h","72h")

level_vector <- c("24h","72h")

e <- ggplot(data_plot[data_plot$time_point_rna_seq %in% level_vector,], aes(x = factor(time_point_rna_seq,level_vector), y = log2FC)) +
  
  geom_violin(
    aes(color = condition), width = 0.3, size = 0.2,
    position = position_dodge(0.8)) + theme_bw() + ylim(-3,1) +
  
  geom_boxplot(
    aes(color = condition), width = 0.05, size = 0.2,
    position = position_dodge(0.8)) + theme_bw() + ylim(-3,1) + theme(legend.position = "none")

ggsave("violin-BRG-BRMdtag-neg-log2FC1-atac-chip-all-rna-seq-supenh.pdf",height = 4, width = 5, dpi=300)
print(e)
dev.off()








### CC1 CC2dtag

### reading data and merging

data_atac <- read.csv("differential_analysis.deseq_result.all_comparisons.csv")

data_sup_enh <- read.table("superenhancer_consensus_overlap.bed",header = F,sep = "\t")

colnames(data_sup_enh) <- c("chr","pos1","pos2")

data_sup_enh$index <- paste(data_sup_enh$chr,":",data_sup_enh$pos1,"-",data_sup_enh$pos2,sep = "")

data_atac_sup <- merge(data_atac,data_sup_enh,"index")


data_want <- data_atac_sup[data_atac_sup$comparison_name %in% c("ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_3h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_6h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_24h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_72h_vs_combined_control"),]

data_want_sig <- data_want[data_want$padj < 0.01 & data_want$log2FoldChange < -1,]

data_chip <- read.table("macs2-log2-FC-rpm-based.txt",header = T)

data_chip_want <- data_chip[,c(1,28,31,34)]

data_atac_chip <- merge(data_want_sig,data_chip_want,"index")

data_gene <- read.csv("BAF_Timecourse_Set3.gene.annotation.singleline.txt")

colnames(data_gene)[1] <- "index"

data_rna <- read.table("rna-seq-all-comp-consensus-log2FConly.txt",header = T)

colnames(data_rna)[2] <- "gene_name"

data_rna_gene <- merge(data_gene,data_rna,"gene_name")

data_rna_gene_want <- data_rna_gene[,c("gene_name","index","gene_id","log2FoldChange_CC1_CC2_F1_dTAG_47_2h_vs_HAP1_CC1_CC2_F1_DMSO","log2FoldChange_CC1_CC2_F1_dTAG_47_3h_vs_HAP1_CC1_CC2_F1_DMSO","log2FoldChange_CC1_CC2_F1_dTAG_47_6h_vs_HAP1_CC1_CC2_F1_DMSO","log2FoldChange_CC1_CC2_F1_dTAG_47_24h_vs_HAP1_CC1_CC2_F1_DMSO","log2FoldChange_CC1_CC2_F1_dTAG_47_72h_vs_HAP1_CC1_CC2_F1_DMSO")]

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

data_plot$time_point_rna_seq[1:256] <- "2h"

data_plot$time_point_rna_seq[257:512] <- "3h"

data_plot$time_point_rna_seq[513:768] <- "6h"

data_plot$time_point_rna_seq[769:1024] <- "24h"

data_plot$time_point_rna_seq[1025:1280] <- "72h"

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"]))))

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"]))))

wilcox.test(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"])),unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"])))

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"]))))

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"]))))

wilcox.test(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"])),unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"])))



### Making violin plots

level_vector <- c("24h","72h")

e <- ggplot(data_plot[data_plot$time_point_rna_seq %in% level_vector,], aes(x = factor(time_point_rna_seq,level_vector), y = log2FC)) +
  
  geom_violin(
    aes(color = condition), width = 0.3, size = 0.2,
    position = position_dodge(0.8)) + theme_bw() + ylim(-3,1) +
  
  geom_boxplot(
    aes(color = condition), width = 0.05, size = 0.2,
    position = position_dodge(0.8)) + theme_bw() + ylim(-3,1) + theme(legend.position = "none")

ggsave("violin-CC1-CC2dtag-neg-log2FC1-atac-chip-all-rna-seq-supenhancer.pdf",height = 4, width = 5, dpi=300)
print(e)
dev.off()