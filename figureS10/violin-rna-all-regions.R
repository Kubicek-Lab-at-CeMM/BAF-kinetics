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

data_rna_gene_want <- data_rna_gene[,c("gene_name","index","gene_id","log2FoldChange_WT_SMARCA4_F1_dTAG_47_2h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_3h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_6h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_24h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_72h_vs_HAP1_WT_SMARCA4_F1_DMSO")]

data_atac_chip_rna <- merge(data_atac_chip,data_rna_gene_want,"index")

### subsetting 

data_atac_chip_rna_down <- data_atac_chip_rna[data_atac_chip_rna$FC_WT_SMARCA4dTAG_DMSO_24h_H3K27ac < -1, ]

data_atac_chip_rna_down$condition <- "dtag-atac-decrease-24h-H3K27ac-chip-decrease"

print(length(data_atac_chip_rna_down$condition))

data_atac_chip_rna_no_response <- data_atac_chip_rna[data_atac_chip_rna$FC_WT_SMARCA4dTAG_DMSO_24h_H3K27ac > -1 & data_atac_chip_rna$FC_WT_SMARCA4dTAG_DMSO_24h_H3K27ac < 1, ]

data_atac_chip_rna_no_response$condition <- "dtag-atac-decrease-24h-H3K27ac-chip-no-response"

print(length(data_atac_chip_rna_no_response$condition))

data_plot_merge <- rbind(data_atac_chip_rna_down,data_atac_chip_rna_no_response)

### rna seq row wise

data_frame_list <- list(data_plot_merge[,c("log2FoldChange_WT_SMARCA4_F1_dTAG_47_2h_vs_HAP1_WT_SMARCA4_F1_DMSO","condition")],data_plot_merge[,c("log2FoldChange_WT_SMARCA4_F1_dTAG_47_3h_vs_HAP1_WT_SMARCA4_F1_DMSO","condition")],data_plot_merge[,c("log2FoldChange_WT_SMARCA4_F1_dTAG_47_6h_vs_HAP1_WT_SMARCA4_F1_DMSO","condition")],data_plot_merge[,c("log2FoldChange_WT_SMARCA4_F1_dTAG_47_24h_vs_HAP1_WT_SMARCA4_F1_DMSO","condition")],data_plot_merge[,c("log2FoldChange_WT_SMARCA4_F1_dTAG_47_72h_vs_HAP1_WT_SMARCA4_F1_DMSO","condition")])

data_plot <- rbindlist(data_frame_list,use.names=F)

colnames(data_plot)[1] <- "log2FC"

data_plot$time_point_rna_seq <- NA

data_plot$time_point_rna_seq[1:2171] <- "2h"

data_plot$time_point_rna_seq[2172:4342] <- "3h"

data_plot$time_point_rna_seq[4343:6513] <- "6h"

data_plot$time_point_rna_seq[6514:8684] <- "24h"

data_plot$time_point_rna_seq[8685:10855] <- "72h"

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
  aes(color = condition), show.legend = FALSE, width = 0.3, size = 0.2,
  position = position_dodge(0.8)) + theme_bw()  + ylim(-6.5,5) +
  
  geom_boxplot(
    aes(color = condition), show.legend = FALSE, width = 0.05, size = 0.2,
    position = position_dodge(0.8)) + theme_bw() + ylim(-6.5,5) + theme(legend.position = "none")

ggsave("violin-WT-SMARCA4dtag-neg-log2FC1-atac-chip-rna-seq-late-tp.pdf",height = 4, width = 5, dpi=300)
print(e)
dev.off()




### BRG BRMdtag

### reading data and merging

data_atac <- read.csv("differential_analysis.deseq_result.all_comparisons.csv")

data_want <- data_atac[data_atac$comparison_name %in% c("ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_72h_vs_DMSO"),]

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

print(length(data_atac_chip_rna_down$condition))

data_atac_chip_rna_no_response <- data_atac_chip_rna[data_atac_chip_rna$FC_SMARCA4KO_SMARCA2dTAG_DMSO_24h_H3K27ac > -1 & data_atac_chip_rna$FC_SMARCA4KO_SMARCA2dTAG_DMSO_24h_H3K27ac < 1, ]

data_atac_chip_rna_no_response$condition <- "dtag-atac-decrease-24h-H3K27ac-chip-no-response"

print(length(data_atac_chip_rna_no_response$condition))

data_plot_merge <- rbind(data_atac_chip_rna_down,data_atac_chip_rna_no_response)

### rna seq row wise

data_frame_list <- list(data_plot_merge[,c("log2FoldChange_BRG_BRM_C8_dTAG_47_2h_vs_HAP1_BRG_BRM_C8_DMSO","condition")],data_plot_merge[,c("log2FoldChange_BRG_BRM_C8_dTAG_47_3h_vs_HAP1_BRG_BRM_C8_DMSO","condition")],data_plot_merge[,c("log2FoldChange_BRG_BRM_C8_dTAG_47_6h_vs_HAP1_BRG_BRM_C8_DMSO","condition")],data_plot_merge[,c("log2FoldChange_BRG_BRM_C8_dTAG_47_24h_vs_HAP1_BRG_BRM_C8_DMSO","condition")],data_plot_merge[,c("log2FoldChange_BRG_BRM_C8_dTAG_47_72h_vs_HAP1_BRG_BRM_C8_DMSO","condition")])

data_plot <- rbindlist(data_frame_list,use.names=F)

colnames(data_plot)[1] <- "log2FC"

data_plot$time_point_rna_seq <- NA

data_plot$time_point_rna_seq[1:986] <- "2h"

data_plot$time_point_rna_seq[987:1972] <- "3h"

data_plot$time_point_rna_seq[1973:2958] <- "6h"

data_plot$time_point_rna_seq[2959:3944] <- "24h"

data_plot$time_point_rna_seq[3945:4930] <- "72h"

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"]))))

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"]))))

wilcox.test(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"])),unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"])))

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"]))))

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"]))))

wilcox.test(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"])),unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"])))



### Making box plots


level_vector <- c("24h","72h")

e <- ggplot(data_plot[data_plot$time_point_rna_seq %in% level_vector,], aes(x = factor(time_point_rna_seq,level_vector), y = log2FC)) +
  
  geom_violin(
    aes(color = condition), show.legend = FALSE, width = 0.3, size = 0.2,
    position = position_dodge(0.8)) + theme_bw()  + ylim(-6.5,5) +
  
  geom_boxplot(
    aes(color = condition), show.legend = FALSE, width = 0.05, size = 0.2,
    position = position_dodge(0.8)) + theme_bw() + ylim(-6.5,5) + theme(legend.position = "none")

ggsave("violin-BRG-BRMdtag-neg-log2FC1-atac-chip-rna-seq-late-tp.pdf",height = 4, width = 5, dpi=300)
print(e)
dev.off()








### CC1 CC2dtag

### reading data and merging

data_atac <- read.csv("/Users/sgrosche/Development/191127-BAF-timecourse-batch3/results/differential_analysis/differential_analysis.deseq_result.all_comparisons.215254.csv")

data_want <- data_atac[data_atac$comparison_name %in% c("ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_3h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_6h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_24h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_72h_vs_combined_control"),]

data_want_sig <- data_want[data_want$padj < 0.01 & data_want$log2FoldChange < -1,]

data_chip <- read.table("/Users/sgrosche/Development/200526-Chip-Seq-Batch1/plots/heatmap_correlation_cov_atac_regions_w_atac_FC_rpm_based/200629-macs2-log2-FC-rpm-based.txt",header = T)

data_chip_want <- data_chip[,c(1,28,31,34)]

data_atac_chip <- merge(data_want_sig,data_chip_want,"index")

data_gene <- read.csv("/Users/sgrosche/Development/191127-BAF-timecourse-batch3/results/191203-BAF_Timecourse_Set3.gene.annotation.singleline.txt")

colnames(data_gene)[1] <- "index"

data_rna <- read.table("/Users/sgrosche/Development/200610-RNA-Seq/data/200803-rna-seq-all-comp-consensus-log2FConly.txt",header = T)

colnames(data_rna)[2] <- "gene_name"

data_rna_gene <- merge(data_gene,data_rna,"gene_name")

#write.table(data_rna_gene,"200818-rna-seq-all-comp-consensus-log2FConly-w-index-dtag-atac.txt",quote = F,sep = "\t",row.names = F)

data_rna_gene_want <- data_rna_gene[,c("gene_name","index","gene_id","log2FoldChange_CC1_CC2_F1_dTAG_47_2h_vs_HAP1_CC1_CC2_F1_DMSO","log2FoldChange_CC1_CC2_F1_dTAG_47_3h_vs_HAP1_CC1_CC2_F1_DMSO","log2FoldChange_CC1_CC2_F1_dTAG_47_6h_vs_HAP1_CC1_CC2_F1_DMSO","log2FoldChange_CC1_CC2_F1_dTAG_47_24h_vs_HAP1_CC1_CC2_F1_DMSO","log2FoldChange_CC1_CC2_F1_dTAG_47_72h_vs_HAP1_CC1_CC2_F1_DMSO")]

data_atac_chip_rna <- merge(data_atac_chip,data_rna_gene_want,"index")

#write.table(data_atac_chip_rna,"200819-CC1-CC2-dtag-atac-neg-log2FC1-w-chip-rna-data.txt",quote = F,sep = "\t",row.names = F)

### subsetting 

data_atac_chip_rna_down <- data_atac_chip_rna[data_atac_chip_rna$FC_SMARCC1KO_SMARCC2dTAG_DMSO_24h_H3K27ac < -1, ]

data_atac_chip_rna_down$condition <- "dtag-atac-decrease-24h-H3K27ac-chip-decrease"

print(length(data_atac_chip_rna_down$condition))

data_atac_chip_rna_no_response <- data_atac_chip_rna[data_atac_chip_rna$FC_SMARCC1KO_SMARCC2dTAG_DMSO_24h_H3K27ac > -1 & data_atac_chip_rna$FC_SMARCC1KO_SMARCC2dTAG_DMSO_24h_H3K27ac < 1, ]

data_atac_chip_rna_no_response$condition <- "dtag-atac-decrease-24h-H3K27ac-chip-no-response"

print(length(data_atac_chip_rna_no_response$condition))

data_plot_merge <- rbind(data_atac_chip_rna_down,data_atac_chip_rna_no_response)

### rna seq row wise

data_frame_list <- list(data_plot_merge[,c("log2FoldChange_CC1_CC2_F1_dTAG_47_2h_vs_HAP1_CC1_CC2_F1_DMSO","condition")],data_plot_merge[,c("log2FoldChange_CC1_CC2_F1_dTAG_47_3h_vs_HAP1_CC1_CC2_F1_DMSO","condition")],data_plot_merge[,c("log2FoldChange_CC1_CC2_F1_dTAG_47_6h_vs_HAP1_CC1_CC2_F1_DMSO","condition")],data_plot_merge[,c("log2FoldChange_CC1_CC2_F1_dTAG_47_24h_vs_HAP1_CC1_CC2_F1_DMSO","condition")],data_plot_merge[,c("log2FoldChange_CC1_CC2_F1_dTAG_47_72h_vs_HAP1_CC1_CC2_F1_DMSO","condition")])

data_plot <- rbindlist(data_frame_list,use.names=F)

colnames(data_plot)[1] <- "log2FC"

data_plot$time_point_rna_seq <- NA

data_plot$time_point_rna_seq[1:3409] <- "2h"

data_plot$time_point_rna_seq[3410:6818] <- "3h"

data_plot$time_point_rna_seq[6819:10227] <- "6h"

data_plot$time_point_rna_seq[10228:13636] <- "24h"

data_plot$time_point_rna_seq[13637:17045] <- "72h"

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"]))))

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"]))))

wilcox.test(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"])),unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "24h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"])))

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"]))))

print(median(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"]))))

wilcox.test(unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-no-response"),"log2FC"])),unlist(na.omit(data_plot[(data_plot$time_point_rna_seq == "72h" & data_plot$condition == "dtag-atac-decrease-24h-H3K27ac-chip-decrease"),"log2FC"])),exact = T)


### Making violin plots

#level_vector <- c("2h","3h","6h","24h","72h")

level_vector <- c("24h","72h")

e <- ggplot(data_plot[data_plot$time_point_rna_seq %in% level_vector,], aes(x = factor(time_point_rna_seq,level_vector), y = log2FC)) +
  
  geom_violin(
    aes(color = condition), show.legend = FALSE, width = 0.3, size = 0.2,
    position = position_dodge(0.8)) + theme_bw() + ylim(-6.5,5) + 
  
  geom_boxplot(
    aes(color = condition), show.legend = FALSE, width = 0.05, size = 0.2,
    position = position_dodge(0.8)) + theme_bw() + ylim(-6.5,5) + theme(legend.position = "none")

ggsave("200929-violin-CC1-CC2dtag-neg-log2FC1-atac-chip-rna-seq-late-tp-same-axis-wo-legend.pdf",height = 4, width = 5, dpi=300)
print(e)
dev.off()
