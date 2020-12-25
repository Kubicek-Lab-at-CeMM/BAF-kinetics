### installing packages

library(ggplot2)

### Setting working directory

setwd("../data/")

### Reading in data

data_diff <- read.table("differential_analysis.deseq_result.all_comparisons.csv",header = T,sep = ",")

for (element in unique(data_diff$comparison_name)){

name = element

data_plot <- data_diff[data_diff$comparison_name == name,]

data_plot$neg_log10_p <- -log10(data_plot$padj)

data_plot$coloring <- ifelse(data_plot$padj < 0.01 & abs(data_plot$log2FoldChange) > 1,'sig',"not sig")


volcano_plot <- ggplot(data_plot) + geom_point(aes(x= log2FoldChange, y=neg_log10_p, color = coloring)) + 
   scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() +
  ylim(0,35) 

pdf(paste("volcano-plot-",name,".pdf"))
print(volcano_plot)
dev.off()

data_sig <- na.omit(data_plot[data_plot$coloring== "sig",])

write.table(data_sig,paste("sig-genes-",element,".txt",sep = ""),quote = F,sep = "\t",row.names = F)

data_sig_up <- data_sig[data_sig[,"log2FoldChange"] > 0,]

write.table(data_sig_up,paste("sig-upregulated-genes-",element,"-atac.txt",sep = ""),quote = F,sep = "\t",row.names = F)

data_sig_down <- data_sig[data_sig[,"log2FoldChange"] < 0,]

write.table(data_sig_down,paste("sig-downregulated-genes-",element,"-atac.txt",sep = ""),quote = F,sep = "\t",row.names = F)



}




sample_vector <- c("ATAC-seq_HAP1_A2_A4_F5_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_A2_A4_F5_dTAG_47_72h_vs_DMSO","ATAC-seq_HAP1_ARID2KO_vs_DMSO","ATAC-seq_HAP1_ARID2KO_vs_WT","ATAC-seq_HAP1_BRGKO_vs_WT")

for (element in sample_vector){
  
  name = element
  
  data_plot <- data_diff[data_diff$comparison_name == name,]
  
  data_plot$neg_log10_p <- -log10(data_plot$padj)
  
  data_plot$coloring <- ifelse(data_plot$padj < 0.01 & abs(data_plot$log2FoldChange) > 1,'sig',"not sig")
  
  
  volcano_plot <- ggplot(data_plot) + geom_point(aes(x= log2FoldChange, y=neg_log10_p, color = coloring)) + 
    scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() +
    ylim(0,20) 
  
  pdf(paste("volcano-plot-",name,"axis-per-cellline.pdf"))
  print(volcano_plot)
  dev.off()
  
}





sample_vector <- c("ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_24h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_2h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_3h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_6h_vs_DMSO","ATAC-seq_HAP1_BRG_BRM_C8_dTAG_47_72h_vs_DMSO","ATAC-seq_HAP1_BRGKO_vs_DMSO","ATAC-seq_HAP1_BRGKO_vs_WT","ATAC-seq_HAP1_BRMKO_vs_WT")

for (element in sample_vector){
  
  name = element
  
  data_plot <- data_diff[data_diff$comparison_name == name,]
  
  data_plot$neg_log10_p <- -log10(data_plot$padj)
  
  data_plot$coloring <- ifelse(data_plot$padj < 0.01 & abs(data_plot$log2FoldChange) > 1,'sig',"not sig")
  
  
  volcano_plot <- ggplot(data_plot) + geom_point(aes(x= log2FoldChange, y=neg_log10_p, color = coloring)) + 
    scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() +
    ylim(0,20) 
  
  pdf(paste("volcano-plot-",name,"axis-per-cellline.pdf"))
  print(volcano_plot)
  dev.off()
  
}




sample_vector <- c("ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_3h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_6h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_24h_vs_combined_control","ATAC-seq_HAP1_CC1_CC2_F1_dTAG_47_72h_vs_combined_control","ATAC-seq_HAP1_CC1KO_vs_combined_control","ATAC-seq_HAP1_CC1KO_vs_WT","ATAC-seq_HAP1_CC2_vs_WT")

for (element in sample_vector){
  
  name = element
  
  data_plot <- data_diff[data_diff$comparison_name == name,]
  
  data_plot$neg_log10_p <- -log10(data_plot$padj)
  
  data_plot$coloring <- ifelse(data_plot$padj < 0.01 & abs(data_plot$log2FoldChange) > 1,'sig',"not sig")
  
  
  volcano_plot <- ggplot(data_plot) + geom_point(aes(x= log2FoldChange, y=neg_log10_p, color = coloring)) + 
    scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() +
    ylim(0,35) 
  
  pdf(paste("volcano-plot-",name,"axis-per-cellline.pdf"))
  print(volcano_plot)
  dev.off()
  
}
