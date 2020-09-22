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

### for the first 5 clusters

### macs2 data

### reading macs2 rpm based FC data and qn count data

data_macs2_qn <- read.csv("Chip.matrix_qn.norm.for.atac.consensus.batch3.regions.csv")

colnames(data_macs2_qn)[1] <- "index"


data_macs2_rpm_FC <- read.table("macs2-log2-FC-rpm-based.txt",header = T)


data_macs2_rpm_FC <- data_macs2_rpm[,c("index","FC_WT_SMARCA4dTAG_DMSO_24h_ARID1A","FC_WT_SMARCA4dTAG_DMSO_24h_BRD4","FC_WT_SMARCA4dTAG_DMSO_24h_H3K27ac","FC_SMARCA4KO_SMARCA2dTAG_DMSO_24h_ARID1A","FC_SMARCA4KO_SMARCA2dTAG_DMSO_24h_BRD4","FC_SMARCA4KO_SMARCA2dTAG_DMSO_24h_H3K27ac","FC_SMARCC1KO_SMARCC2dTAG_DMSO_24h_ARID1A","FC_SMARCC1KO_SMARCC2dTAG_DMSO_24h_BRD4","FC_SMARCC1KO_SMARCC2dTAG_DMSO_24h_H3K27ac")]

data_macs2_complete <- merge(data_macs2_qn,data_macs2_rpm_FC,"index")

write.table(data_macs2_complete,"macs2-FC-rpm-based-qn-counts-chip-batch1-merged.txt",quote = F,sep = "\t",row.names = F)

columns_wanted <- c("index","WT_SMARCA4dTAG_DMSO_ARID1A","WT_SMARCA4dTAG_dTAG47_24h_ARID1A","WT_SMARCA4dTAG_DMSO_BRD4","WT_SMARCA4dTAG_dTAG47_24h_BRD4","WT_SMARCA4dTAG_DMSO_H3K27ac","WT_SMARCA4dTAG_dTAG47_24h_H3K27ac","FC_WT_SMARCA4dTAG_DMSO_24h_ARID1A","FC_WT_SMARCA4dTAG_DMSO_24h_BRD4","FC_WT_SMARCA4dTAG_DMSO_24h_H3K27ac")

data_macs2_plot <- data_macs2_complete[,columns_wanted]

rowsd_macs2_1<-apply(data_macs2_plot[,c(2:7)],1,sd)
rowmean_macs2_1<-rowMeans(data_macs2_plot[,c(2:7)])
Ztrans_macs2_1<-(data_macs2_plot[,c(2:7)]-rowmean_macs2_1)/rowsd_macs2_1

Ztrans_macs2_1$index <- data_macs2_plot$index

FC_macs2 <- data_macs2_plot[,c(1,8:10)]


### reading in cluster data FC based

cluster1 <- read.table("Cluster1-timecourse-only-regions-3-clusters-FC-based.txt",header = T)

cluster1$index <- row.names(cluster1)

cluster2 <- read.table("Cluster2-timecourse-only-regions-3-clusters-FC-based.txt",header = T)

cluster2$index <- row.names(cluster2)

cluster3 <- read.table("Cluster3-timecourse-only-regions-3-clusters-FC-based.txt",header = T)

cluster3$index <- row.names(cluster3)

data_merge <- rbind(cluster1,cluster2,cluster3)

### reading atac data

data_norm <- fread("BAF_Timecourse.matrix_cqn_norm_merged_replicates.tsv",data.table = F)

SMARCA4vector <- c("WT_DMSO","WT_SMARCA4_F1_DMSO","WT_SMARCA4_F1_dTAG_47_2h","WT_SMARCA4_F1_dTAG_47_3h","WT_SMARCA4_F1_dTAG_47_6h","WT_SMARCA4_F1_dTAG_47_24h","WT_SMARCA4_F1_dTAG_47_72h","SMARCA4_KO_DMSO")

data_norm_SMARCA4 <- data_norm[,c("index",SMARCA4vector)]

rowsd<-apply(data_norm_SMARCA4[,-1],1,sd)
rowmean<-rowMeans(data_norm_SMARCA4[,-1])
Ztrans<-(data_norm_SMARCA4[,-1]-rowmean)/rowsd

Ztrans$index <- data_norm_SMARCA4$index

merge_all <- merge(data_merge,Ztrans,"index")

### merging data

merge_all_chip <- merge(merge_all,Ztrans_macs2_1,"index")

merge_all_chip_2 <- merge(merge_all_chip,FC_macs2,"index")

row.names(merge_all_chip_2) <- merge_all_chip_2$index

colnames(merge_all_chip_2) <- gsub("X","",colnames(merge_all_chip_2))

colnames(merge_all_chip_2)[2:7] <- paste("FC_",colnames(merge_all_chip_2)[2:7],sep="")




### plotting heatmap

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

# with own hclust function
d <- dist(merge_all_chip_2[,-c(1,7:24)], method = "canberra") # distance matrix
fit <- hclust(d, method="ward.D")

## only timecourse 3 clusters

ht1 <- Heatmap(as.matrix(merge_all_chip_2[,2:6]), name = "log2 FC or Z-score cqn counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),split = 3,cluster_rows = fit,cluster_columns = F,col=col_fun)

row_order_want <- unlist(row_order(ht1))

ht2 <- Heatmap(as.matrix(merge_all_chip_2[,7]),column_title = "", row_title = "",column_labels = "FC_KO_SMARCA4",column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_columns = FALSE,cluster_rows = FALSE,col=col_fun,show_heatmap_legend = F,row_order = row_order_want)

ht3 <- Heatmap(as.matrix(merge_all_chip_2[,8:15]),column_title = "", row_title = "", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_columns = FALSE,cluster_rows = FALSE,col=col_fun,show_heatmap_legend = F,row_order = row_order_want)

ht4 <- Heatmap(as.matrix(merge_all_chip_2[,16:21]),column_title = "", row_title = "", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_columns = FALSE,cluster_rows = FALSE,col=col_fun,show_heatmap_legend = F,row_order = row_order_want)

ht5 <- Heatmap(as.matrix(merge_all_chip_2[,22:24]),column_title = "", row_title = "", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_columns = FALSE,cluster_rows = FALSE,col=col_fun,show_heatmap_legend = F,row_order = row_order_want)


pdf("Heatmap-FC-clustered-canberra-ward-rows-SMARCA4-dtag-tc.pdf")
draw(ht1+ht3+ht1+ht2+ht4+ht5)
dev.off()
