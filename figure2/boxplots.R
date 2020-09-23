### loading packages

library(ggplot2)
library(ggsci)
library(MultNonParam)
if(!require(BSDA)){install.packages("BSDA")}
if(!require(DescTools)){install.packages("DescTools")}
library(BSDA)
library(lawstat)

### setting working directory

setwd("../data/")

### reading in data

data <- read.table("FC_p_allclusters_SMARCA4_timecourse_KO_differential_regions.txt",header = T)

### writing N regions

write.table(table(data$cluster)/6,"N-regions-clusters.txt",quote = F,sep = "\t",row.names = F,col.names = T)

### making boxplots

cluster_vector <- c("cluster1_cluster2_KO","cluster3_KO","cluster1_tc","cluster2_tc","cluster3_tc")

for (element in cluster_vector){
  
  # Basic box plot
  
  cluster_name = element
  
  data_filtered <- data[data$cluster == paste(cluster_name,sep = ""),]
  
  level_vector <- c("WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO","WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO","WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO","WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO","WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO","BRGKO_vs_WT")
  
  
  p <- ggplot(data_filtered, aes(x=factor(comparison_name,level_vector), y=log2FoldChange, color = comparison_name)) + 
    geom_boxplot()  + ylim(-9,4) + theme_bw() + theme(legend.position = "none",axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_colour_npg() + geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)
  
  
  ggsave(paste("Boxplot-accessibility-log2FC-SMARCA4-timecourse-KO-",cluster_name,"-wolegend.pdf",sep = ""), plot = p, device = "pdf", path = NULL,
         scale = 1, height= 6, width = 8, dpi = 300, limitsize = TRUE)
  
  
}
