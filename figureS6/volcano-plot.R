### installing packages

library(corrplot)
library(data.table)
library(factoextra)
library(NbClust)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)

### setting working directory

setwd("../data/")

### vulcano plots 

data_diff <- read.csv2("differential_analysis.deseq_result.all_comparisons_novartis_inhibitor.csv",stringsAsFactors = F,sep = ",")

data_diff[2:(length(colnames(data_diff))-1)] <- apply(data_diff[2:(length(colnames(data_diff))-1)],2,function(x) as.numeric(x))

data_diff_sig <- data_diff[data_diff$padj < 0.01 & abs(data_diff$log2FoldChange) > 1,]

write.table(data_diff_sig,"sig-diff-regions-only-all-comp.txt",quote = F,sep = "\t",row.names = F)

for (element in unique(data_diff$comparison_name)){
  
  name = element
  print(name)
  
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
