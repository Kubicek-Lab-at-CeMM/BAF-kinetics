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
library(ggsci)

### setting working directory

setwd("../data/")

### for the SL clusters

### macs2 data

### reading macs2 qn and rpm data

data_macs2_qn <- read.csv("Chip.matrix_qn.norm.for.atac.consensus.batch3.regions.csv")

colnames(data_macs2_qn)[1] <- "index"

data_macs2_rpm_FC <- read.table("macs2-log2-FC-rpm-based.txt",header = T)

### making Z macs2 qn

rowsd<-apply(data_macs2_qn[,-c(1,2,6,10,14,18,22)],1,sd)
rowmean<-rowMeans(data_macs2_qn[,-c(1,2,6,10,14,18,22)])
Z_macs2_qn<-(data_macs2_qn[,-c(1,2,6,10,14,18,22)]-rowmean)/rowsd

Z_macs2_qn$index <- data_macs2_qn$index


### chosing samples to plot

data_macs2_rpm_FC <- data_macs2_rpm_FC[,c("index","FC_WT_SMARCA4dTAG_DMSO_24h_ARID1A","FC_WT_SMARCA4dTAG_DMSO_24h_BRD4","FC_WT_SMARCA4dTAG_DMSO_24h_H3K27ac","FC_SMARCA4KO_SMARCA2dTAG_DMSO_24h_ARID1A","FC_SMARCA4KO_SMARCA2dTAG_DMSO_24h_BRD4","FC_SMARCA4KO_SMARCA2dTAG_DMSO_24h_H3K27ac","FC_SMARCC1KO_SMARCC2dTAG_DMSO_24h_ARID1A","FC_SMARCC1KO_SMARCC2dTAG_DMSO_24h_BRD4","FC_SMARCC1KO_SMARCC2dTAG_DMSO_24h_H3K27ac")]

data_macs2_complete <- merge(data_macs2_qn,data_macs2_rpm_FC,"index")

columns_wanted <- c("index","FC_WT_SMARCA4dTAG_DMSO_24h_ARID1A","FC_WT_SMARCA4dTAG_DMSO_24h_BRD4","FC_WT_SMARCA4dTAG_DMSO_24h_H3K27ac","FC_SMARCA4KO_SMARCA2dTAG_DMSO_24h_ARID1A","FC_SMARCA4KO_SMARCA2dTAG_DMSO_24h_BRD4","FC_SMARCA4KO_SMARCA2dTAG_DMSO_24h_H3K27ac","FC_SMARCC1KO_SMARCC2dTAG_DMSO_24h_ARID1A","FC_SMARCC1KO_SMARCC2dTAG_DMSO_24h_BRD4","FC_SMARCC1KO_SMARCC2dTAG_DMSO_24h_H3K27ac")

data_macs2_plot <- data_macs2_complete[,colnames(data_macs2_complete) %in% columns_wanted]


### reading in Sl clusters

data_merge <- read.table("cluster_SL.annotated.bed")

colnames(data_merge) <- c("index","cluster")

### Reading in cqn counts

data_norm <- fread("BAF_Timecourse.matrix_cqn_norm_merged_replicates.tsv",data.table = F)

data_norm_want <- data_norm[,c(1,24:30,16,8,2:7,9,10:15,23,22,17:21)]

rowsd<-apply(data_norm_want[,-1],1,sd)
rowmean<-rowMeans(data_norm_want[,-1])
Ztrans<-(data_norm_want[,-1]-rowmean)/rowsd

Ztrans$index <- data_norm_want$index

merge_all <- merge(data_merge,Ztrans,"index")

### reading fold change data

data_FC <- read.csv("differential_analysis.deseq_result.all_comparisons.csv")

### getting FC table

FC_table_list <- list()

comparison_vector <- c("ATAC-seq_HAP1_WT_vs_SMARCA4_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO","ATAC-seq_HAP1_BRGKO_vs_WT",
                       "ATAC-seq_HAP1_ARID2KO_vs_WT","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_72h_vs_DMSO","ATAC-seq_HAP1_BRMKO_vs_WT",
                       "ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_72h_vs_DMSO",
                       "ATAC-seq_HAP1_CC2_vs_WT","ATAC-seq_HAP1_CC1KO_vs_WT","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_3h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_6h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_24h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_72h_vs_combined_control")

for (comparison in comparison_vector){
  
  data_interest <- data_FC[data_FC$comparison == paste(comparison,sep=""),c("index","log2FoldChange")]
  
  colnames(data_interest)[2] <- paste(colnames(data_interest)[2],"_",comparison,sep = "")
  
  FC_table_list[[comparison]] <- data_interest
  
}

FC_table <- Reduce(function(x,y) merge(x,y,"index"),FC_table_list)

FC_table$`log2FoldChange_ATAC-seq_HAP1_SMARCA4_DMSO_vs_WT` <- -(FC_table$`log2FoldChange_ATAC-seq_HAP1_WT_vs_SMARCA4_DMSO`)

write.table(FC_table,"consensus-FC-atac-data-per-column.txt",quote = F,sep = "\t",row.names = F)

### merging data

merge_all_fc <- merge(merge_all,FC_table,"index")

merge_all_chip <- merge(merge_all_fc,data_macs2_plot,"index")

merge_all_chip_2 <- merge(merge_all_chip,Z_macs2_qn,"index")

row.names(merge_all_chip_2) <- merge_all_chip_2$index

merge_all_chip_2$index <- NULL

merge_all_chip_2$cluster <- NULL

### plotting heatmap

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

# with own hclust function
d <- dist(merge_all_chip_2[,-c(1:37,43,49:50,55:82)], method = "canberra") # distance matrix
fit <- hclust(d, method="ward.D")

### fold change heatmap

ht1 <- Heatmap(as.matrix(merge_all_chip_2[,c(31:35,36,37:54)]), name = "log2 FC", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),split = 6,cluster_rows = fit,cluster_columns = F,col=col_fun)

row_order_want <- unlist(row_order(ht1))

ht2 <- Heatmap(as.matrix(merge_all_chip_2[,56:64]),column_title = "", row_title = "",column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_columns = FALSE,cluster_rows = FALSE,col=col_fun,show_heatmap_legend = F,row_order = row_order_want)


pdf("Heatmap-FC-clustered-SL-tc.pdf")
draw(ht1+ht2)
dev.off()
