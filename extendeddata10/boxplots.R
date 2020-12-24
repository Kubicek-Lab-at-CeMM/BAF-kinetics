### setting working directory

setwd("../data/")

### reading in FC data

data <- read.table("200302-supenhancer-diff-any-condition-merged-deseq2-output.txt",header = T)

### making boxplots

### SMARCA4dtag 

level_vector <- c("ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO","ATAC-seq_HAP1_BRGKO_vs_WT")

data_filtered <- data[data$comparison_name %in% level_vector,]

p <- ggplot(data_filtered, aes(x=factor(comparison_name,level_vector), y=log2FoldChange, color = comparison_name)) + 
  geom_boxplot() + ylim(-10,5) + theme_bw() + theme(axis.text.x = element_blank(),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_colour_npg() + geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)

ggsave(paste("boxplot-accessibility-log2FC-SMARCA4-timecourse-KO-diff-supenh-wolegend.pdf",sep = ""), plot = p, device = "pdf", path = NULL,
       scale = 1, height= 6, width = 8, dpi = 300, limitsize = TRUE)

q <- ggplot(data_filtered, aes(x=factor(comparison_name,level_vector), y=log2FoldChange, color = comparison_name)) + 
  geom_boxplot() + ylim(-10,5) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_colour_npg() + geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)

ggsave(paste("boxplot-accessibility-log2FC-SMARCA4-timecourse-KO-diff-supenh-withlegend.pdf",sep = ""), plot = q, device = "pdf", path = NULL,
       scale = 1, height= 6, width = 8, dpi = 300, limitsize = TRUE)



### ARID2KO SMARCA4dtag 

level_vector <- c("ATAC-seq_HAP1_ARID2KO_vs_WT","ATAC-seq_HAP1_BRGKO_vs_WT","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_72h_vs_DMSO")

data_filtered <- data[data$comparison_name %in% level_vector,]

p <- ggplot(data_filtered, aes(x=factor(comparison_name,level_vector), y=log2FoldChange, color = comparison_name)) + 
  geom_boxplot() + ylim(-10,5) + theme_bw() + theme(axis.text.x = element_blank(),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_colour_npg() + geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)

ggsave(paste("boxplot-accessibility-log2FC-ARID2KO-SMARCA4dtag-timecourse-KO-diff-supenh-wolegend.pdf",sep = ""), plot = p, device = "pdf", path = NULL,
       scale = 1, height= 6, width = 8, dpi = 300, limitsize = TRUE)


q <- ggplot(data_filtered, aes(x=factor(comparison_name,level_vector), y=log2FoldChange, color = comparison_name)) + 
  geom_boxplot() + ylim(-10,5) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_colour_npg() + geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)

ggsave(paste("boxplot-accessibility-log2FC-ARID2KO-SMARCA4dtag-timecourse-KO-diff-supenh-withlegend.pdf",sep = ""), plot = q, device = "pdf", path = NULL,
       scale = 1, height= 6, width = 8, dpi = 300, limitsize = TRUE)



### CC1 CC2dtag 

level_vector <- c("ATAC-seq_HAP1_CC1KO_vs_WT","ATAC-seq_HAP1_CC2_vs_WT","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_3h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_6h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_24h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_72h_vs_combined_control")

data_filtered <- data[data$comparison_name %in% level_vector,]

p <- ggplot(data_filtered, aes(x=factor(comparison_name,level_vector), y=log2FoldChange, color = comparison_name)) + 
  geom_boxplot() + ylim(-10,5) + theme_bw() + theme(axis.text.x = element_blank(),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_colour_npg() + geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)

ggsave(paste("boxplot-accessibility-log2FC-CC1KO-CC2-timecourse-KO-diff-supenh-wolegend.pdf",sep = ""), plot = p, device = "pdf", path = NULL,
       scale = 1, height= 6, width = 8, dpi = 300, limitsize = TRUE)

q <- ggplot(data_filtered, aes(x=factor(comparison_name,level_vector), y=log2FoldChange, color = comparison_name)) + 
  geom_boxplot() + ylim(-10,5) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_colour_npg() + geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)

ggsave(paste("boxplot-accessibility-log2FC-CC1KO-CC2-timecourse-KO-diff-supenh-withlegend.pdf",sep = ""), plot = q, device = "pdf", path = NULL,
       scale = 1, height= 6, width = 8, dpi = 300, limitsize = TRUE)



### BRG BRMdtag 

level_vector <- c("ATAC-seq_HAP1_BRGKO_vs_WT","ATAC-seq_HAP1_BRMKO_vs_WT","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_72h_vs_DMSO")

data_filtered <- data[data$comparison_name %in% level_vector,]

p <- ggplot(data_filtered, aes(x=factor(comparison_name,level_vector), y=log2FoldChange, color = comparison_name)) + 
  geom_boxplot() + ylim(-10,5) + theme_bw() + theme(axis.text.x = element_blank(),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_colour_npg() + geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)

ggsave(paste("boxplot-accessibility-log2FC-BRGKO-BRM-timecourse-KO-diff-supenh-wolegend.pdf",sep = ""), plot = p, device = "pdf", path = NULL,
       scale = 1, height= 6, width = 8, dpi = 300, limitsize = TRUE)

q <- ggplot(data_filtered, aes(x=factor(comparison_name,level_vector), y=log2FoldChange, color = comparison_name)) + 
  geom_boxplot() + ylim(-10,5) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_colour_npg() + geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)

ggsave(paste("boxplot-accessibility-log2FC-BRGKO-BRM-timecourse-KO-diff-supenh-withlegend.pdf",sep = ""), plot = q, device = "pdf", path = NULL,
       scale = 1, height= 6, width = 8, dpi = 300, limitsize = TRUE)
