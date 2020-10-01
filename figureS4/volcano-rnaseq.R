
### installing packages

library(ggplot2)
library(ggrepel)

### Setting working directory

setwd("../data/")

### Reading in data

data <- read.table("rna-seq-all-comp-consensus.txt",header = T)

sample_vector <- c("WT_SMARCA4_F1_dTAG_47_24h_vs_HAP1_WT_SMARCA4_F1_DMSO","WT_SMARCA4_F1_dTAG_47_2h_vs_HAP1_WT_SMARCA4_F1_DMSO","WT_SMARCA4_F1_dTAG_47_3h_vs_HAP1_WT_SMARCA4_F1_DMSO","WT_SMARCA4_F1_dTAG_47_6h_vs_HAP1_WT_SMARCA4_F1_DMSO","WT_SMARCA4_F1_dTAG_47_72h_vs_HAP1_WT_SMARCA4_F1_DMSO")

for (element in sample_vector){
  
  data_want <- data[, grepl(element,colnames(data)) ]
    
  colnames(data_want)[2] <- "gene_name"
  
  data_want$neg_log10_p <- -log10(data_want[,paste("pvalue_",element,sep = "")])
  
  data_want$coloring <- ifelse(data_want[,paste("padj_",element,sep = "")] < 0.01 & abs(data_want[,paste("log2FoldChange_",element,sep = "")]) > 1,'sig',"not sig")
  
  data_want$name <- ifelse(data_want$coloring == "sig",as.character(data_want$gene_name),"")
  
  volcano_plot <- ggplot(data_want) + geom_point(aes(x= data_want[,paste("log2FoldChange_",element,sep = "")], y=neg_log10_p, color = coloring)) + 
    scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() +
      xlim(-4,4)
    
  
  pdf(paste("volcano-plot-",element,"-deseq2-results.pdf",sep=""))
  print(volcano_plot)  
  dev.off()
  
}


element = "WT_DMSO_vs_HAP1_SMARCA4_KO_DMSO"

data_want <- data[, grepl(element,colnames(data)) ]

colnames(data_want)[2] <- "gene_name"

data_want$neg_log10_p <- -log10(data_want[,paste("pvalue_",element,sep = "")])

data_want$neg_log10_p_clean <- ifelse(data_want$neg_log10_p == -Inf,NA,data_want$neg_log10_p)

data_want$neg_log10_p_clean <- ifelse(data_want$neg_log10_p == Inf,NA,data_want$neg_log10_p)

data_want$coloring <- ifelse(data_want[,paste("padj_",element,sep = "")] < 0.01 & abs(data_want[,paste("log2FoldChange_",element,sep = "")]) > 1,'sig',"not sig")

data_want$name <- ifelse(data_want$coloring == "sig",as.character(data_want$gene_name),"")

volcano_plot <- ggplot(data_want) + geom_point(aes(x= data_want[,paste("log2FoldChange_",element,sep = "")], y=neg_log10_p_clean, color = coloring)) + 
  scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() + xlim(-10,10)


pdf(paste("volcano-plot-",element,"-deseq2-results.pdf",sep=""))
print(volcano_plot)  
dev.off()
