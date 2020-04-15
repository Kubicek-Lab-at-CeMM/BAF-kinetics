### installing packages

library(ggplot2)

### Setting working directory

setwd("../data/")

### Reading in data

data_diff <- read.table("differential_analysis.deseq_result.all_comparisons.csv",header = T,sep = ",")

sample_vector <- c("ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO","ATAC-seq_HAP1_BRGKO_vs_WT")

for (element in sample_vector){

name = element

data_plot <- data_diff[data_diff$comparison_name == name,]

data_plot$neg_log10_p <- -log10(data_plot$padj)

data_plot$coloring <- ifelse(data_plot$padj < 0.01 & abs(data_plot$log2FoldChange) > 1,'sig',"not sig")


volcano_plot <- ggplot(data_plot) + geom_point(aes(x= log2FoldChange, y=neg_log10_p, color = coloring)) + 
   scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() +
  ylim(0,20) + xlim(-8,8)

pdf(paste("volcano-plot-",name,".pdf"))
print(volcano_plot)
dev.off()

}
