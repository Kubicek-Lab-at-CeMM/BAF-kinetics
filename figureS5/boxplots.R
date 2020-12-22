library(ggsci)
library(ggplot2)

### setting working directory

setwd("../data/")

### reading in data

data_FC <- read.table("ATP-depletion-timecourse-allregions-new-deseq2.txt",header = T,stringsAsFactors = F)

colnames(data_FC)[1] <- "index_compound"

### reading merge file

merge_file <- read.table("dtag_ATP_overlap_merge.bed",header=F,stringsAsFactors = F)

colnames(merge_file) <- c("chr_dtag","pos1_dtag","pos2_dtag","chr_compound","pos1_compound","pos2_compound")

merge_file$index_dtag <- paste(merge_file$chr_dtag,":",merge_file$pos1_dtag,"-",merge_file$pos2_dtag,sep="")

merge_file$index_compound <- paste(merge_file$chr_compound,":",merge_file$pos1_compound,"-",merge_file$pos2_compound,sep="")

### merging data

data_FC_merge <- merge(data_FC,merge_file,"index_compound")

### making boxplots per old cluster

level_vector <- c("5min_Novartis_inhibitor","10min_Novartis_inhibitor","30min_Novartis_inhibitor","1h_Novartis_inhibitor","6h_Novartis_inhibitor","24h_Novartis_inhibitor","6h_BI_vs_WT_DMSO","24h_BI_vs_WT_DMSO")

cluster_vector <- c("cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8","cluster9","cluster10","cluster11")

for (cluster in cluster_vector){
 
  cluster_data <- read.table(paste("",cluster,"-annotated.bed",sep=""),header = F)
  
  colnames(cluster_data) <- c("index_dtag","cluster")
  
  data_FC_merge_cluster <- merge(cluster_data,data_FC_merge,"index_dtag")
  
  print(table(data_FC_merge_cluster$comparison_name))
  
  p <- ggplot(data_FC_merge_cluster, aes(x=factor(comparison_name,level_vector), y=log2FoldChange, color = comparison_name)) + 
    geom_boxplot() + ylim(-9,6) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_colour_npg() + geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)
  
  ggsave(paste("boxplot-accessibility-log2FC-",cluster,".pdf",sep = ""), plot = p, device = "pdf", path = NULL,
         scale = 1, height= 6, width = 8, dpi = 300, limitsize = TRUE)
  
     
}


