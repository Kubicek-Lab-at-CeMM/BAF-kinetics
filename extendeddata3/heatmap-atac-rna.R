### installing packages

library(data.table)
library(factoextra)
library(NbClust)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(ggplot2)
library(grid)
library(gridExtra)

### setting working directory

setwd("../data/")

### reading data 

data_atac <- read.table("diff-expression-data-atac-seq-tc-column-based.txt",header = T)

data_gene <- read.table("rna-seq-all-comp-consensus-log2FConly.txt",header = T)

colnames(data_gene)[2] <- "gene_name"

data_gene$log2FoldChange_HAP1_SMARCA4_KO_DMSO_vs_WT_DMSO <- -data_gene$log2FoldChange_WT_DMSO_vs_HAP1_SMARCA4_KO_DMSO

data_gene$log2FoldChange_HAP1_SMARCC1_KO_DMSO_vs_WT_DMSO <- -data_gene$log2FoldChange_WT_DMSO_vs_HAP1_SMARCC1_KO_DMSO

data_gene$log2FoldChange_HAP1_SMARCC2_KO_DMSO_vs_WT_DMSO <- -data_gene$log2FoldChange_WT_DMSO_vs_HAP1_SMARCC2_KO_DMSO

gene_anno <- read.csv("BAF_Timecourse_Set3.gene.annotation.singleline.txt",header = T)

cluster_regions <- read.table("cluster_all.annotated.bed")

colnames(cluster_regions) <- c("index","cluster")

### merging data

data_merge_gene_anno <- merge(data_gene,gene_anno,"gene_name")

colnames(data_merge_gene_anno)[46] <- "index"

data_merge_atac_gene <- merge(data_atac,data_merge_gene_anno,"index")

data_merge_cluster <- merge(cluster_regions,data_merge_atac_gene,"index")


### WT SMARCA4dtag time-course

cluster_vector <- c("cluster1","cluster2","cluster3")

ht_list <- list()

for (cluster in cluster_vector){

print(cluster)

data_merge_cluster_want <- data_merge_cluster[data_merge_cluster$cluster == cluster,]

data_plot <- as.matrix(data_merge_cluster_want[,c("log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO_atac","log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO_atac","log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO_atac","log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO_atac","log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO_atac","log2FoldChange_ATAC.seq_HAP1_BRGKO_vs_WT_atac", "log2FoldChange_WT_SMARCA4_F1_dTAG_47_2h_vs_HAP1_WT_SMARCA4_F1_DMSO", "log2FoldChange_WT_SMARCA4_F1_dTAG_47_3h_vs_HAP1_WT_SMARCA4_F1_DMSO", "log2FoldChange_WT_SMARCA4_F1_dTAG_47_6h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_24h_vs_HAP1_WT_SMARCA4_F1_DMSO", "log2FoldChange_WT_SMARCA4_F1_dTAG_47_72h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_HAP1_SMARCA4_KO_DMSO_vs_WT_DMSO")])

print(length(data_plot[,1]))

row.names(data_plot) <- paste(data_merge_cluster_want$index,"_",row.names(data_plot),sep = "")

### making heatmap

col_fun_1 = colorRamp2(c(-1.5, 0, 1.5), c("blue", "gray90", "red"))


## with own hclust function
dist_1 <- dist(data_plot[,c(7:11)], method = "manhattan") # distance matrix
fit_1 <- hclust(dist_1, method= "ward.D")


ht1 <- Heatmap(data_plot[,c(7:11)], name = "log2 FC", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = fit_1,cluster_columns = F,col=col_fun_1)

row_ordering <- unlist(row_order(ht1))

ht2 <- Heatmap(data_plot[,c(1:5)], name = "log2 FC", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = F,cluster_columns = F,col=col_fun_1,row_order = row_ordering)

ht3 <- Heatmap(data_plot[,c(6)], column_labels = "log2FoldChange_ATAC.seq_HAP1_BRGKO_vs_WT_atac", name = "log2 FC", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = F,cluster_columns = F,col=col_fun_1,row_order = row_ordering)

ht4 <- Heatmap(data_plot[,c(7:11)], name = "log2 FC", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = F,cluster_columns = F,col=col_fun_1,row_order = row_ordering)

ht5 <- Heatmap(data_plot[,c(12)], column_labels = "log2FoldChange_HAP1_SMARCA4_KO_DMSO_vs_WT_DMSO" ,name = "log2 FC", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = F,cluster_columns = F,col=col_fun_1,row_order = row_ordering)

plot <- ht1 + ht2 + ht3 + ht4 + ht5

pdf(paste("heatmap-atac-WT-SMARCA4dtag-",cluster,"clusters-w-rna-seq.pdf",sep=""))
draw(ht1 + ht2 + ht3 + ht4 + ht5)
dev.off()

ht_list[[cluster]] <- plot

}



### WT SMARCA4 KO

cluster_vector <- c("cluster4","cluster5")

for (cluster in cluster_vector){
  
  print(cluster)
  
  data_merge_cluster_want <- data_merge_cluster[data_merge_cluster$cluster == cluster,]
  
  data_plot <- as.matrix(data_merge_cluster_want[,c("log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO_atac","log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO_atac","log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO_atac","log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO_atac","log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO_atac","log2FoldChange_ATAC.seq_HAP1_BRGKO_vs_WT_atac", "log2FoldChange_WT_SMARCA4_F1_dTAG_47_2h_vs_HAP1_WT_SMARCA4_F1_DMSO", "log2FoldChange_WT_SMARCA4_F1_dTAG_47_3h_vs_HAP1_WT_SMARCA4_F1_DMSO", "log2FoldChange_WT_SMARCA4_F1_dTAG_47_6h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_24h_vs_HAP1_WT_SMARCA4_F1_DMSO", "log2FoldChange_WT_SMARCA4_F1_dTAG_47_72h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_HAP1_SMARCA4_KO_DMSO_vs_WT_DMSO")])
  
  print(length(data_plot[,1]))
  
  row.names(data_plot) <- paste(data_merge_cluster_want$index,"_",row.names(data_plot),sep = "")
  
  ### making heatmap
  
  col_fun_1 = colorRamp2(c(-1.5, 0, 1.5), c("blue", "gray90", "red"))
  
  
  ## with own hclust function
  dist_1 <- dist(data_plot[,c(7:12)], method = "manhattan") # distance matrix
  fit_1 <- hclust(dist_1, method= "ward.D")
  
  
  ht1 <- Heatmap(data_plot[,c(7:11)], name = "log2 FC", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = fit_1,cluster_columns = F,col=col_fun_1)
  
  row_ordering <- unlist(row_order(ht1))
  
  ht2 <- Heatmap(data_plot[,c(1:5)], name = "log2 FC", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = F,cluster_columns = F,col=col_fun_1,row_order = row_ordering)
  
  ht3 <- Heatmap(data_plot[,c(6)], column_labels = "log2FoldChange_ATAC.seq_HAP1_BRGKO_vs_WT_atac", name = "log2 FC", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = F,cluster_columns = F,col=col_fun_1,row_order = row_ordering)
  
  ht4 <- Heatmap(data_plot[,c(7:11)], name = "log2 FC", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = F,cluster_columns = F,col=col_fun_1,row_order = row_ordering)
  
  ht5 <- Heatmap(data_plot[,c(12)], column_labels = "log2FoldChange_HAP1_SMARCA4_KO_DMSO_vs_WT_DMSO" ,name = "log2 FC", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = F,cluster_columns = F,col=col_fun_1,row_order = row_ordering)
  
  plot <- ht1 + ht2 + ht3 + ht4 + ht5
  
  pdf(paste("heatmap-atac-WT-SMARCA4KO-",cluster,"clusters-w-rna-seq-nosplit.pdf",sep=""))
  draw(ht1 + ht2 + ht3 + ht4 + ht5)
  dev.off()
  
  ht_list[[cluster]] <- plot
  
}

