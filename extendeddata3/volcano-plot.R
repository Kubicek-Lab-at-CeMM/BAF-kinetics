### installing packages

library(ggplot2)
library(ggrepel)

### Setting working directory

setwd("../data/")

### Reading in data

sample_vector <- c("10min-Novartis-vs-3h-WT-DMSO","1h-Novartis-vs-3h-WT-DMSO","30min-Novartis-vs-3h-WT-DMSO","3h-dtag-WT-SMARCA4-dtag-vs-3h-DMSO-WT-SMARCA4-dtag","3h-dtag-BRG_BRM-dtag-vs-3h-DMSO-BRG_BRM-dtag")

for (element in sample_vector){

  data <- read.table(paste("deseq2-",element,"-gene-full-genebody.txt",sep=""),header = T)
  
  data$neg_log10_p <- -log10(data$pvalue)
  
  data$coloring <- ifelse(data$padj < 0.01 & abs(data$log2FoldChange) > 1,'sig',"not sig")
  
  data$name <- ifelse(data$coloring == "sig",as.character(data$gene),"")
  
  volcano_plot <- ggplot(data) + geom_point(aes(x= log2FoldChange, y=neg_log10_p, color = coloring)) + 
    scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() +
    ylim(0,15) + xlim(-6,6) + geom_text_repel(aes(x= log2FoldChange, y=neg_log10_p,label=name),size = 3)
  
  pdf(paste("volcano-plot-gene-level-full-body-",element,"-deseq2-results.pdf",sep=""))
  print(volcano_plot)  
  dev.off()
  
  data_clean <- na.omit(data)
  
  data_sig <- data_clean[data_clean$coloring == "sig",]
  
  write.table(data_sig,paste("data-sig-deseq2-gene-level-full-body-",element,".txt",sep=""),quote = F, sep = "\t",row.names = F)
  
}


