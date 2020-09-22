### getting packages

library(data.table)
library(ggplot2)

### setting working directory

setwd("../data/")

### functions

### calculating median and max/min per sample

median_maker <- function(x) {
  
  clean <- na.omit(x)
  
  wo_inf <- clean[is.infinite(clean)==F]
  
  median(wo_inf)
}

se_maker <- function(x) {
  
  clean <- na.omit(x)
  
  wo_inf <- clean[is.infinite(clean)==F]
  
  sd(wo_inf) / sqrt(length(wo_inf)) 

}


### reading in foldchange data

### atac ATP depletion time-course + pro seq

FC_atac_ATP <- read.table("ATP-depletion-data-new-deseq2-treatment-per-column-all-regions.txt",header = T)

FC_atac_ATP_tc <- FC_atac_ATP[,c("index_compound","log2FoldChange_N_5min","log2FoldChange_N_10min","log2FoldChange_N_30min","log2FoldChange_N_1h","log2FoldChange_N_6h","log2FoldChange_N_24h")]

FC_pro_seq_ATP <- fread("deseq2-Novartis-proseq-merge-column.txt",data.table = F)

FC_pro_seq_ATP_want <- FC_pro_seq_ATP[,c("index","log2FoldChange_10min_Novartis","log2FoldChange_30min_Novartis","log2FoldChange_1h_Novartis")]

colnames(FC_pro_seq_ATP_want)[1:4] <- c("index_compound",paste(colnames(FC_pro_seq_ATP_want)[2:4],"_proseq",sep=""))

merge_ATP_depl_atac_proseq <- merge(FC_atac_ATP_tc,FC_pro_seq_ATP_want,"index_compound")

### merging fold change data

data_merge_2 <- merge_ATP_depl_atac_proseq

### looping over clusters

cluster_vector <- c("clusterI","clusterII","clusterIII","clusterIV","clusterV")

for (element in cluster_vector){
  
  
  cluster = read.table(paste("",element,"-FC-clustered.bed",sep=""),header = F)
  
  colnames(cluster) <- c("chr","pos1","pos2")
  
  cluster$index_compound <- paste(cluster$chr,":",cluster$pos1,"-",cluster$pos2,sep = "")
  
  data_merge_cluster_2 <- merge(data_merge_2,cluster,"index_compound")

  ### calculating median and max/min per sample

  median_2 <- apply(data_merge_cluster_2[2:10],2,median_maker)
  
  se_2_pos <- median_2 + apply(data_merge_cluster_2[2:10],2,se_maker)
  
  se_2_neg <- median_2 - (apply(data_merge_cluster_2[2:10],2,se_maker))
  
  median_all <- median_2
  
  se_pos_all <- se_2_pos
  
  se_neg_all <- se_2_neg
  
  
  data <- data.frame(median_all,se_pos_all,se_neg_all)
  
  data$type <- NA
  
  data$timepoint <- NA
  
  data$type[1:6] <- "atac-Novartis-inhibitor"
  
  data$type[7:9] <- "proseq-Novartis-inhibitor"
  
  data$timepoint[1:9] <- c("5min","10min","30min","1h","6h","24h","10min","30min","1h")
  
  
  # Make the plot
  
  data$timepoint <- factor(data$timepoint, levels=c("5min","10min","30min","1h","6h","24h"))
  
  pdf(paste("lineplot-",element,".pdf",sep=""))
  print(ggplot(data=data, aes(x=timepoint, y=median_all, group=type)) +
          geom_line(aes(color=type)) +
          geom_point(aes(color=type), size=2) +
          geom_ribbon(aes(ymin = se_neg_all, ymax = se_pos_all, fill = type),alpha = 0.3) +
          ylab("Fold change") + theme_bw() + ylim(-2.0,1.5) + theme(axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
  dev.off()
  
}








