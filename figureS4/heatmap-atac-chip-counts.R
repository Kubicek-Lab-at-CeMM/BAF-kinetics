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

### reading macs2 qn count data chip batch 1

data_macs2_qn_batch1 <- read.csv("Chip.matrix_qn.norm.for.atac.consensus.batch3.regions.csv")

colnames(data_macs2_qn_batch1)[1] <- "index"

columns_wanted <- c("index","WT_SMARCA4dTAG_DMSO_ARID1A","WT_SMARCA4dTAG_dTAG47_24h_ARID1A","WT_SMARCA4dTAG_DMSO_BRD4","WT_SMARCA4dTAG_dTAG47_24h_BRD4","WT_SMARCA4dTAG_DMSO_H3K27ac","WT_SMARCA4dTAG_dTAG47_24h_H3K27ac")

data_macs2_plot <- data_macs2_qn_batch1[,columns_wanted]

rowsd_macs2_1<-apply(data_macs2_plot[,c(2:7)],1,sd)
rowmean_macs2_1<-rowMeans(data_macs2_plot[,c(2:7)])
Ztrans_macs2_1<-(data_macs2_plot[,c(2:7)]-rowmean_macs2_1)/rowsd_macs2_1

Ztrans_macs2_1$index <- data_macs2_plot$index

colnames(Ztrans_macs2_1)[1:6] <- paste(colnames(Ztrans_macs2_1)[1:6],"_batch1",sep = "")

### reading macs2 qn count data chip batch 2

data_macs2_qn_batch2 <- read.csv("Chip-Set2-quantile-norm-atac-dtag-consensus.csv")

colnames(data_macs2_qn_batch2)[1] <- "index"

columns_wanted <- c("index","WT_SMARCA4dTAG_DMSO_H3K27ac","WT_SMARCA4dTAG_3h_dTAG_H3K27ac","WT_SMARCA4dTAG_6h_dTAG_H3K27ac","WT_SMARCA4dTAG_24h_dTAG_H3K27ac","WT_SMARCA4dTAG_72h_dTAG_H3K27ac","SMARCA4KO_DMSO_H3K27ac")

data_macs2_plot_batch2 <- data_macs2_qn_batch2[,columns_wanted]

rowsd_macs2_2<-apply(data_macs2_plot_batch2[,c(2:7)],1,sd)
rowmean_macs2_2<-rowMeans(data_macs2_plot_batch2[,c(2:7)])
Ztrans_macs2_2<-(data_macs2_plot_batch2[,c(2:7)]-rowmean_macs2_2)/rowsd_macs2_2

Ztrans_macs2_2$index <- data_macs2_plot_batch2$index

colnames(Ztrans_macs2_2)[1:6] <- paste(colnames(Ztrans_macs2_2)[1:6],"_batch2",sep = "")


Ztrans_overall <- merge(Ztrans_macs2_1,Ztrans_macs2_2,"index")


### reading atac data

comparison_vector <- c("ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO","ATAC-seq_HAP1_BRGKO_vs_WT")

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


### reading in cluster data FC based

cluster1 <- read.table("cluster1-timecourse-only-regions-3-clusters-FC-based.txt",header = T)

cluster1$index <- row.names(cluster1)

cluster2 <- read.table("cluster2-timecourse-only-regions-3-clusters-FC-based.txt",header = T)

cluster2$index <- row.names(cluster2)

cluster3 <- read.table("cluster3-timecourse-only-regions-3-clusters-FC-based.txt",header = T)

cluster3$index <- row.names(cluster3)

data_merge <- rbind(cluster1,cluster2,cluster3)

cluster4 <- read.table("cluster1-SMARCA4KO-only-regions-3-clusters-FC-based.txt",header = T)

cluster4$index <- row.names(cluster4)

cluster5 <- read.table("cluster2-SMARCA4KO-only-regions-3-clusters-FC-based.txt",header = T)

cluster5$index <- row.names(cluster5)

cluster6 <- read.table("cluster3-SMARCA4KO-only-regions-3-clusters-FC-based.txt",header = T)

cluster6$index <- row.names(cluster6)

data_merge_KO <- rbind(cluster4,cluster5,cluster6)

### getting data to plot dtag tc clusters

data_atac_chip <- merge(atac_data_clean,Ztrans_overall,"index")

data_atac_chip_3cluster <- merge(data_atac_chip,data_merge,"index")

data_atac_chip_3cluster_plot <- as.matrix(data_atac_chip_3cluster[,c(2:19)])

row.names(data_atac_chip_3cluster_plot) <- data_atac_chip_3cluster$index

### plotting heatmap

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

# with own hclust function
d <- dist(data_atac_chip_3cluster_plot[,c(1:5)], method = "canberra") # distance matrix
fit <- hclust(d, method="ward.D")

## only timecourse 3 clusters

ht1 <- Heatmap(data_atac_chip_3cluster_plot[,1:5], name = "log2 FC or Z-score qn counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),split = 3,cluster_rows = fit,cluster_columns = F,col=col_fun)

row_order_want <- unlist(row_order(ht1))

ht2 <- Heatmap(data_atac_chip_3cluster_plot[,6],column_title = "", row_title = "",column_labels = "log2FoldChange_ATAC-seq_HAP1_BRGKO_vs_WT",column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_columns = FALSE,cluster_rows = FALSE,col=col_fun,show_heatmap_legend = F,row_order = row_order_want)

ht3 <- Heatmap(data_atac_chip_3cluster_plot[,7:12],column_title = "", row_title = "", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_columns = FALSE,cluster_rows = FALSE,col=col_fun,show_heatmap_legend = F,row_order = row_order_want)

ht4 <- Heatmap(data_atac_chip_3cluster_plot[,13:18],column_title = "", row_title = "", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_columns = FALSE,cluster_rows = FALSE,col=col_fun,show_heatmap_legend = F,row_order = row_order_want)

cluster1 <- data_atac_chip_3cluster_plot[row_order(ht1)[[1]],]

cluster1_correct_order <- data.frame(row.names(cluster1))

colnames(cluster1_correct_order) <- "index"

cluster2 <- data_atac_chip_3cluster_plot[row_order(ht1)[[2]],]

cluster2_correct_order <- data.frame(row.names(cluster2))

colnames(cluster2_correct_order) <- "index"

cluster3 <- data_atac_chip_3cluster_plot[row_order(ht1)[[3]],]

cluster3_correct_order <- data.frame(row.names(cluster3))

colnames(cluster3_correct_order) <- "index"

cluster_all <- rbind(cluster1_correct_order,cluster2_correct_order,cluster3_correct_order)

pdf("Heatmap-FC-clustered-SMARCA4-dtag-tc-chip-qn-counts-Zscore.pdf")
draw(ht1+ht2+ht3+ht4)
dev.off()


### getting data to plot KO clusters

data_atac_chip_3cluster <- merge(data_atac_chip,data_merge_KO,"index")

data_atac_chip_3cluster_plot <- as.matrix(data_atac_chip_3cluster[,c(2:19)])

row.names(data_atac_chip_3cluster_plot) <- data_atac_chip_3cluster$index

### plotting heatmap

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

# with own hclust function
d <- dist(data_atac_chip_3cluster_plot[,c(1:6)], method = "canberra") # distance matrix
fit <- hclust(d, method="ward.D")

## only timecourse 3 clusters

ht1 <- Heatmap(data_atac_chip_3cluster_plot[,1:5], name = "log2 FC or Z-score qn counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),split = 3,cluster_rows = fit,cluster_columns = F,col=col_fun)

row_order_want <- unlist(row_order(ht1))

ht2 <- Heatmap(data_atac_chip_3cluster_plot[,6],column_title = "", row_title = "",column_labels = "log2FoldChange_ATAC-seq_HAP1_BRGKO_vs_WT",column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_columns = FALSE,cluster_rows = FALSE,col=col_fun,show_heatmap_legend = F,row_order = row_order_want)

ht3 <- Heatmap(data_atac_chip_3cluster_plot[,7:12],column_title = "", row_title = "", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_columns = FALSE,cluster_rows = FALSE,col=col_fun,show_heatmap_legend = F,row_order = row_order_want)

ht4 <- Heatmap(data_atac_chip_3cluster_plot[,13:18],column_title = "", row_title = "", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_columns = FALSE,cluster_rows = FALSE,col=col_fun,show_heatmap_legend = F,row_order = row_order_want)

cluster1 <- data_atac_chip_3cluster_plot[row_order(ht1)[[1]],]

cluster1_correct_order <- data.frame(row.names(cluster1))

colnames(cluster1_correct_order) <- "index"

cluster2 <- data_atac_chip_3cluster_plot[row_order(ht1)[[2]],]

cluster2_correct_order <- data.frame(row.names(cluster2))

colnames(cluster2_correct_order) <- "index"

cluster3 <- data_atac_chip_3cluster_plot[row_order(ht1)[[3]],]

cluster3_correct_order <- data.frame(row.names(cluster3))

colnames(cluster3_correct_order) <- "index"

cluster_all <- rbind(cluster1_correct_order,cluster2_correct_order,cluster3_correct_order)

pdf("Heatmap-FC-clustered-SMARCA4-KO-chip-qn-counts-Zscore.pdf")
draw(ht1+ht2+ht3+ht4)
dev.off()
