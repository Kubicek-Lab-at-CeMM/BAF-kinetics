### installing packages
pkgs <- c("factoextra",  "NbClust")
install.packages(pkgs)
library(data.table)
library(factoextra)
library(NbClust)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)

### setting working directory

setwd("../data/")

### reading in FC based clusters for SMARCA4KO

cluster1 <- read.table("Cluster1-SMARCA4KO-only-regions-3-clusters-FC-based.txt",header = T)

cluster1$index <- row.names(cluster1)

cluster2 <- read.table("Cluster2-SMARCA4KO-only-regions-3-clusters-FC-based.txt",header = T)

cluster2$index <- row.names(cluster2)

cluster3 <- read.table("Cluster3-SMARCA4KO-only-regions-3-clusters-FC-based.txt",header = T)

cluster3$index <- row.names(cluster3)

### merging data

data_merge <- rbind(cluster1,cluster2,cluster3)

### Reading in cqn counts

data_norm <- fread("BAF_Timecourse.matrix_cqn_norm_merged_replicates.tsv",data.table = F)

SMARCA4vector <- c("WT_DMSO","WT_SMARCA4_F1_DMSO","WT_SMARCA4_F1_dTAG_47_2h","WT_SMARCA4_F1_dTAG_47_3h","WT_SMARCA4_F1_dTAG_47_6h","WT_SMARCA4_F1_dTAG_47_24h","WT_SMARCA4_F1_dTAG_47_72h","SMARCA4_KO_DMSO")

data_norm_SMARCA4 <- data_norm[,c("index",SMARCA4vector)]

rowsd<-apply(data_norm_SMARCA4[,-1],1,sd)
rowmean<-rowMeans(data_norm_SMARCA4[,-1])
Ztrans<-(data_norm_SMARCA4[,-1]-rowmean)/rowsd

Ztrans$index <- data_norm_SMARCA4$index

merge_all <- merge(data_merge,Ztrans,"index")

row.names(merge_all) <- merge_all$index

### plotting heatmap

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

## clustering

d <- dist(merge_all[,-c(1,8:15)], method = "canberra") # distance matrix
fit <- hclust(d, method="ward.D")

## only KO 3 clusters

ht1 <- Heatmap(as.matrix(merge_all[,2:6]), name = "log2 FC or cqn counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),split = 3,cluster_rows = fit,cluster_columns = F,col=col_fun)

ht2 <- Heatmap(as.matrix(merge_all[,7]),column_title = "", row_title = "",column_labels = "KO_SMARCA4",column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_columns = FALSE,cluster_rows = FALSE,col=col_fun,show_heatmap_legend = F)

ht3 <- Heatmap(as.matrix(merge_all[,8:15]),column_title = "", row_title = "", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_columns = FALSE,cluster_rows = FALSE,col=col_fun,show_heatmap_legend = F)


pdf("Heatmap-FC-clustered-canberra-ward-rows-SMARCA4-KO-clustered-all-split3-FC-based-w-cqn-counts-next-split.pdf")
draw(ht1+ht3+ht1+ht2)
dev.off()
