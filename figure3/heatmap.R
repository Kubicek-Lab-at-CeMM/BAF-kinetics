### installing packages

library(data.table)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(RColorBrewer)

### setting working directory

setwd("../data/")

### reading in novel diff regions

ARID2KO <- read.table("SL-ARID2KO-bg-SMARCA4-dtag-additional-differential-regions.bed",header = T)

ARID2KO$index <- paste(ARID2KO$chr,":",ARID2KO$pos1,"-",ARID2KO$pos2,sep = "")

BRGKO <- read.table("SL-BRG-bg-BRM-dtag-additional-differential-regions.bed",header=T)

BRGKO$index <- paste(BRGKO$chr,":",BRGKO$pos1,"-",BRGKO$pos2,sep = "")

CC1KO <- read.table("SL-CC1-bg-CC2-dtag-additional-differential-regions.bed",header = T)

CC1KO$index <- paste(CC1KO$chr,":",CC1KO$pos1,"-",CC1KO$pos2,sep = "")

diff_all <- rbind(ARID2KO,BRGKO,CC1KO)

cqn_counts <- read.table("BAF_Timecourse.matrix_cqn_norm_merged_replicates.tsv",header = T)

### merging data

diff_all_unique <- unique(data.frame(diff_all$index))

colnames(diff_all_unique) <- "index"

merge_diff <- merge(diff_all_unique,cqn_counts,"index")

### making heatmap matrix

counts_diff_cqn_sorted <- merge_diff[,c(24:30,16,8,2:7,9,10:15,23,22,17:21)]

counts_diff_cqn_matrix <- as.matrix(counts_diff_cqn_sorted)

row.names(counts_diff_cqn_matrix) <- merge_diff$index

### calculating Zscore

rowsd<-apply(counts_diff_cqn_matrix,1,sd)
rowmean<-rowMeans(counts_diff_cqn_matrix)
Ztrans<-(counts_diff_cqn_matrix-rowmean)/rowsd

### setting color range for heatmap

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

## clustering

d <- dist(Ztrans, method = "canberra") # distance matrix
fit <- hclust(d, method="ward.D")

## plotting heatmap without clusters

ht1 <- Heatmap(Ztrans, name = "Z score normalized counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = fit,cluster_columns = F,col=col_fun)

png("Heatmap-Zscore-cqn-clustered-canberra-ward-rows-all-diff-regions-SL-additional.png",height = 10000, width = 10000, res = 1200)
print(ht1)
dev.off()


## plotting heatmap with clusters

ht2 <- Heatmap(Ztrans, name = "Z score normalized counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), split = 5, column_names_gp = gpar(fontsize = 6),cluster_rows = fit,cluster_columns = F,col=col_fun)

pdf("Heatmap-Zscore-cqn-clustered-canberra-ward-rows-all-diff-regions-SL-additional-split5.pdf")
print(ht2)
dev.off()
