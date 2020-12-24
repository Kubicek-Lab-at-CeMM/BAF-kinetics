### setting working directory

setwd("../data/")

### reading in data

counts_cqn <- read.csv("BAF_Timecourse_Set3.matrix_cqn_norm.csv",header = T)

deseq_results <- read.table("differential_analysis.deseq_result.all_comparisons.csv",header=T,sep=",")

### choosing only differential regions

deseq_diff <- deseq_results[deseq_results$padj < 0.01 & abs(deseq_results$log2FoldChange) > 1, ]

regions_interest <- unique(deseq_diff$index)

counts_diff_cqn <- counts_cqn[counts_cqn$X %in% regions_interest,]

### making heatmap matrix

counts_diff_cqn_sorted <- counts_diff_cqn[,c(59,60,47,53,49,55,50,56,51,57,48,54,52,58,30,31)]

counts_diff_cqn_matrix <- as.matrix(counts_diff_cqn_sorted)

row.names(counts_diff_cqn_matrix) <- counts_diff_cqn$X

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


pdf("Heatmap-Zscore-cqn-clustered-canberra-ward-rows-all-diff-regions-WT-clone.pdf")
print(ht1)
dev.off()

order_want <- row_order(ht1)

### plotting ARID2KO SMARCA4dtag next to WT SMARCA4dtag

### making heatmap matrix

counts_diff_cqn_sorted <- counts_diff_cqn[,c(59,60,47,53,49,55,50,56,51,57,48,54,52,58,30,31,14,15,2,8,4,10,5,11,6,12,3,9,7,13)]

counts_diff_cqn_matrix <- as.matrix(counts_diff_cqn_sorted)

row.names(counts_diff_cqn_matrix) <- counts_diff_cqn$X

### calculating Zscore

rowsd<-apply(counts_diff_cqn_matrix,1,sd)
rowmean<-rowMeans(counts_diff_cqn_matrix)
Ztrans<-(counts_diff_cqn_matrix-rowmean)/rowsd

ht2 <- Heatmap(Ztrans, name = "Z score normalized counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = F,cluster_columns = F,col=col_fun,row_order = order_want)

pdf("Heatmap-Zscore-cqn-clustered-canberra-ward-rows-all-diff-regions-WT-clone-clusterd-ARID2KO-clone-plotted.pdf")
print(ht2)
dev.off()


### seperate heatmaps

d <- dist(Ztrans[,1:16], method = "canberra") # distance matrix
fit <- hclust(d, method="ward.D")


ht3 <- Heatmap(Ztrans[,1:16], name = "Z score normalized counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = fit,cluster_columns = F,col=col_fun)

ht4 <- Heatmap(Ztrans[,c(1,2,17:30,15,16)], name = "Z score normalized counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = F,cluster_columns = F,col=col_fun)


pdf("Heatmap-Zscore-cqn-clustered-canberra-ward-rows-all-diff-regions-WT-clone-clusterd-ARID2KO-clone-plotted-next-wKOWT.pdf")
print(ht3 + ht4)
dev.off()


### seperate ARID2KO clone added

### making heatmap matrix

counts_diff_cqn_sorted <- counts_diff_cqn[,c(14,15,2,8,4,10,5,11,6,12,3,9,7,13)]

counts_diff_cqn_matrix <- as.matrix(counts_diff_cqn_sorted)

row.names(counts_diff_cqn_matrix) <- counts_diff_cqn$X

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

ht5 <- Heatmap(Ztrans, name = "Z score normalized counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = fit,cluster_columns = F,col=col_fun)


pdf("Heatmap-Zscore-cqn-clustered-canberra-ward-rows-all-diff-regions-ARID2KO-clone.pdf")
print(ht5)
dev.off()












### heatmap withouth ARID2KO (added 14/01/20)

order_want <- row_order(ht1)

ht2 <- Heatmap(Ztrans[,c(1:8,17:22,16,24:29,23)], name = "Z score normalized counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = F,cluster_columns = F,col=col_fun,row_order = order_want)

png("200114-Heatmap-Zscore-cqn-clustered-canberra-ward-rows-all-diff-regions-woARID2KO.png",height = 10000, width = 10000, res = 1200)
print(ht2)
dev.off()

pdf("200114-Heatmap-Zscore-cqn-clustered-canberra-ward-rows-all-diff-regions-woARID2KO.pdf")
print(ht2)
dev.off()
