### getting packages

library(data.table)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(RColorBrewer)

### setting working directory

setwd("/Users/sgrosche/Development/191127-BAF-timecourse-batch3/results/differential_analysis/191210-all-timecourse/")

### reading in data

counts_cqn <- read.table("/Users/sgrosche/Development/191127-BAF-timecourse-batch3/results/191127_BAF_Timecourse_Set3.matrix_cqn_norm_merged_replicates_191126_220139.tsv",header = T)

deseq_results <- read.table("/Users/sgrosche/Development/191127-BAF-timecourse-batch3/results/differential_analysis/differential_analysis.deseq_result.all_comparisons.215254.csv",header=T,sep=",")

### choosing only differential regions

deseq_diff <- deseq_results[deseq_results$padj < 0.01 & abs(deseq_results$log2FoldChange) > 1, ]

regions_interest <- unique(deseq_diff$index)

counts_diff_cqn <- counts_cqn[counts_cqn$index %in% regions_interest,]

### making heatmap matrix

counts_diff_cqn_sorted <- counts_diff_cqn[,c(24:30,16,8,2:7,9,10:15,23,22,17:21)]

counts_diff_cqn_matrix <- as.matrix(counts_diff_cqn_sorted)

row.names(counts_diff_cqn_matrix) <- counts_diff_cqn$index

### calculating Zscore

rowsd<-apply(counts_diff_cqn_matrix,1,sd)
rowmean<-rowMeans(counts_diff_cqn_matrix)
Ztrans<-(counts_diff_cqn_matrix-rowmean)/rowsd

### setting color range for heatmap

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

### calculating distance

d <- dist(Ztrans, method = "canberra") # distance matrix
fit <- hclust(d, method="ward.D")

## plotting heatmap

ht1 <- Heatmap(Ztrans, name = "Z score normalized counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = fit,cluster_columns = F,col=col_fun)

png("191210-Heatmap-Zscore-cqn-clustered-canberra-ward-rows-all-diff-regions.png",height = 10000, width = 10000, res = 1200)
print(ht1)
dev.off()

pdf("200131-Heatmap-Zscore-cqn-clustered-canberra-ward-rows-all-diff-regions.pdf")
print(ht1)
dev.off()
