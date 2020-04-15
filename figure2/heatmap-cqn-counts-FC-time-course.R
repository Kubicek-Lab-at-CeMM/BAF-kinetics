### Getting packages

library(data.table)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(RColorBrewer)

### Setting working directory

setwd("../data/")

### Reading in clusters FC based

cluster1 <- read.table("Cluster1-timecourse-only-regions-3-clusters-FC-based.txt",header = T)

cluster1$index <- row.names(cluster1)

cluster2 <- read.table("Cluster2-timecourse-only-regions-3-clusters-FC-based.txt",header = T)

cluster2$index <- row.names(cluster2)

cluster3 <- read.table("Cluster3-timecourse-only-regions-3-clusters-FC-based.txt",header = T)

cluster3$index <- row.names(cluster3)

### Reading in cqn counts

data_norm <- fread("BAF_Timecourse.matrix_cqn_norm_merged_replicates.tsv",data.table = F)

SMARCA4vector <- c("WT_DMSO","WT_SMARCA4_F1_DMSO","WT_SMARCA4_F1_dTAG_47_2h","WT_SMARCA4_F1_dTAG_47_3h","WT_SMARCA4_F1_dTAG_47_6h","WT_SMARCA4_F1_dTAG_47_24h","WT_SMARCA4_F1_dTAG_47_72h","SMARCA4_KO_DMSO")

data_norm_SMARCA4 <- data_norm[,c("index",SMARCA4vector)]

rowsd<-apply(data_norm_SMARCA4[,-1],1,sd)
rowmean<-rowMeans(data_norm_SMARCA4[,-1])
Ztrans<-(data_norm_SMARCA4[,-1]-rowmean)/rowsd

Ztrans$index <- data_norm_SMARCA4$index

### merging data

clustermerge1 <- merge(cluster1,Ztrans,"index")

clustermerge2 <- merge(cluster2,Ztrans,"index")

clustermerge3 <- merge(cluster3,Ztrans,"index")

clustermergeall <- rbind(clustermerge1,clustermerge2,clustermerge3)

row.names(clustermergeall) <- clustermergeall$index

clustermergeall$index <- NULL

clustermergeallmatrix <- as.matrix(clustermergeall)

### plotting heatmap

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

## clustering

d <- dist(clustermergeallmatrix[,-c(6:14)], method = "canberra") # distance matrix
fit <- hclust(d, method="ward.D")

## only timecourse 3 clusters

ht1 <- Heatmap(clustermergeallmatrix[,1:5], name = "log2 FC or cqn counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),split = 3,cluster_rows = fit,cluster_columns = F,col=col_fun)

ht2 <- Heatmap(clustermergeallmatrix[,6],column_title = "", row_title = "",column_labels = "KO_SMARCA4",column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_columns = FALSE,cluster_rows = FALSE,col=col_fun,show_heatmap_legend = F)

ht3 <- Heatmap(clustermergeallmatrix[,7:14],column_title = "", row_title = "", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_columns = FALSE,cluster_rows = FALSE,col=col_fun,show_heatmap_legend = F)


png("Heatmap-FC-clustered-canberra-ward-rows-SMARCA4-timecourse-clustered-all-split3-FC-based-w-cqn-counts-next-split.png",height = 10000, width = 10000, res = 1200)
print(ht1+ht2+ht3)
dev.off()

pdf("Heatmap-FC-clustered-canberra-ward-rows-SMARCA4-timecourse-clustered-all-split3-FC-based-w-cqn-counts-next-split.pdf")
print(ht1+ht3+ht1+ht2)
dev.off()

