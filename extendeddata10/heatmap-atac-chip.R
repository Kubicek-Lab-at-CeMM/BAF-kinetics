### getting packages

library(data.table)
library(factoextra)
library(NbClust)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(plyr)
library(ggplot2)
library(hexbin)
library(viridis)
library(scales)
library(reshape)
library(ggsci)

### setting working directory

setwd("../data/")

### reading atac data

comparison_vector <- c("ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO","ATAC-seq_HAP1_BRGKO_vs_WT","ATAC-seq_HAP1_ARID2KO_vs_WT","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_72h_vs_DMSO","ATAC-seq_HAP1_BRMKO_vs_WT","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_72h_vs_DMSO","ATAC-seq_HAP1_CC2_vs_WT","ATAC-seq_HAP1_CC1KO_vs_WT","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_3h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_6h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_24h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_72h_vs_combined_control")

data <- read.csv("differential_analysis.deseq_result.all_comparisons.csv")

data_atac_want <- data[data$comparison_name %in% comparison_vector,]

list_atac <- list()

for (comparison in comparison_vector) {
  
  print(comparison)
  data_comp <- data_atac_want[data_atac_want$comparison_name == comparison,]
  data_comp_need <- data_comp[,c("index","log2FoldChange")]
  colnames(data_comp_need)[2] <- paste(colnames(data_comp_need)[2],"_",comparison,sep = "")
  list_atac[[comparison]] <- data_comp_need
}

atac_data_clean <- Reduce(function(x,y) merge(x,y,"index"),list_atac)

write.table(atac_data_clean,"atac-log2FC-column-based-all-consensus.txt",quote = F,sep = "\t",row.names = F)

### reading in chip data

data_chip <- read.table("log2FC-rpm-based-chip-seq-batch2-dtag-atac-consensus-regions.txt",header = T)

data_chip_want <- data_chip[,c("index","log2FCWT_SMARCA4dTAG_3h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_6h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_24h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_72h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA2KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_3h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_6h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_24h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_72h_H3K27ac_vs_DMSO")]

### merging data

data_atac_chip <- merge(atac_data_clean,data_chip_want,"index")

write.table(data_atac_chip,"atac-chip-batch2-log2FC-column-based-all-consensus.txt",quote = F,sep = "\t",row.names = F)

data_atac_chip_want <- data_atac_chip[,-c(8:13,20:25)]

### merging with clusterdata

data_cluster <- read.table("all_11_clusters_right_order.txt",header = F)

colnames(data_cluster) <- "index"

data_plot <- join(data_cluster,data_atac_chip_want,"index")

### plotting heatmap without NA


ht1 <- Heatmap(data_plot[,2:13], name = "log2 FC", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 4),cluster_rows = F,cluster_columns = F,col=col_fun)

ht2 <- Heatmap(data_plot[,14:23], name = "log2 FC", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 4),cluster_rows = F,cluster_columns = F,col=col_fun)


pdf("all-cluster-atac-chip-batch2-WTSMARCA4dtag-BRGKOBRMdtag.pdf")
print(ht1 + ht2)
dev.off()

