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

### Setting working directory

setwd("../data/")

### reading differential expression data and cqn normalized count data

data_diff <- fread("differential_analysis.deseq_result.all_comparisons.csv",data.table = F) 

data_norm <- fread("BAF_Timecourse.matrix_cqn_norm_merged_replicates.tsv",data.table = F)

### identifying differential regions based on the SMARCA4 knockdown timecourse only

SMARCA4vector <- c("ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO")

data_diff_SMARCA4 <- data_diff[data_diff$comparison_name %in% SMARCA4vector,]

data_diff_SMARCA4_filtered <- data_diff_SMARCA4[data_diff_SMARCA4$padj < 0.01 & abs(data_diff_SMARCA4$log2FoldChange) > 1,]

SMARCA4_timecourse_diff_vector <- unique(data_diff_SMARCA4_filtered$index)

### merging normalized countdata and SMARCA4 knockdown timecourse differential regions

data_norm_counts_SMRACA4_timecourse <- data_norm[data_norm$index %in% SMARCA4_timecourse_diff_vector,]

### extracting samples for SMARCA4 timecourse heatmap

samplevector <- c("index","WT_DMSO","WT_SMARCA4_F1_DMSO","WT_SMARCA4_F1_dTAG_47_2h","WT_SMARCA4_F1_dTAG_47_3h","WT_SMARCA4_F1_dTAG_47_6h","WT_SMARCA4_F1_dTAG_47_24h","WT_SMARCA4_F1_dTAG_47_72h","SMARCA4_KO_DMSO")

data_norm_counts_SMARCA4_timecourse_samples <- data_norm_counts_SMRACA4_timecourse[,samplevector]

data_norm_counts_SMARCA4_timecourse_samples_matrix <- as.matrix(data_norm_counts_SMARCA4_timecourse_samples[,-1])

row.names(data_norm_counts_SMARCA4_timecourse_samples_matrix) <- data_norm_counts_SMARCA4_timecourse_samples[,1]

### PCA analysis

res.pca <- prcomp(t(data_norm_counts_SMARCA4_timecourse_samples_matrix), scale = TRUE)


# PCA plot

pdf("190225-pca-plot-PC1and2-merged-samples-colored-sample.pdf",width = 7, height = 5,)
fviz_pca_ind(res.pca,
             col.ind.sup = c("#5ac18e", "#0069b4", "#b51865","#6f1a16", "#00aee7", "#762a02", "#fcd116", "#518561"),
             labelsize = 4,
             pointshape = 20,
             col.ind  = c("WT_DMSO","WT_SMARCA4_F1_DMSO","WT_SMARCA4_F1_dTAG_47_2h","WT_SMARCA4_F1_dTAG_47_3h","WT_SMARCA4_F1_dTAG_47_6h","WT_SMARCA4_F1_dTAG_47_24h","WT_SMARCA4_F1_dTAG_47_72h","SMARCA4_KO_DMSO"),
             repel = T,
             label = F
)
dev.off()

