library(ggsci)
library(ggplot2)

### setting working directory

setwd("../data/")

### reading in fold change data

data_FC <- read.table("differential_analysis.deseq_result.all_comparisons_novartis_inhibitor.csv",header = T,sep = ",",stringsAsFactors = F)

### boxplots per cluster

level_vector <- c("5min_Novartis_inhibitor","10min_Novartis_inhibitor","30min_Novartis_inhibitor","1h_Novartis_inhibitor","6h_Novartis_inhibitor","24h_Novartis_inhibitor","6h_BI_protac","24h_BI_protac")

clustervector <- c("clusterI","clusterII","clusterIII","clusterIV","clusterV")

for (element in clustervector){
  
  cluster <- read.table(paste("",element,"-FC-clustered.bed",sep=""),sep = "\t")
  
  colnames(cluster) <- c("chr","pos1","pos2")
  
  cluster$index <- paste(cluster$chr,":",cluster$pos1,"-",cluster$pos2,sep="")

  data_FC_cluster <-  merge(cluster,data_FC,"index")
  
  q <- ggplot(data_FC_cluster, aes(x=factor(comparison_name,level_vector), y=log2FoldChange, color = comparison_name)) + 
    geom_boxplot() + ylim(-7.5,3) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_colour_npg() + geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)
  
  ggsave(paste("boxplot-accessibility-log2FC-BRM014-inhibitor-timecourse-",element,".pdf",sep = ""), plot = q, device = "pdf", path = NULL,
         scale = 1, height= 6, width = 8, dpi = 300, limitsize = TRUE)
 
  }
  