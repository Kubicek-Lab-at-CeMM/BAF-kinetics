### Loading packages

library(data.table)
library(ggsci)

### Setting working directory

setwd("../data/")

### Reading H3K27ac data

data <- fread("ChIP-seq_HAP1_WT_H3K27ac_921-7.hg38.raw.counts.txt",header = F)

colnames(data) <- c("chr","pos1","pos2","counts")

data$index <- paste(data$chr,":",data$pos1,"-",data$pos2,sep = "")

### Analysis for H3K27ac all regions

### Opening list

cluster_list <- list()

### Iterating over clusters

cluster_vector <- c("cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8","cluster9","cluster10")

for (element in cluster_vector){

print(element)
  
cluster <- read.table(paste("/data/",element,".bed",sep = ""),header = F)

colnames(cluster) <- c("chr","pos1","pos2")

cluster$index <- paste(cluster$chr,":",cluster$pos1,"-",cluster$pos2,sep = "")

data_cluster <- merge(cluster,data,"index")

data_cluster$cluster <- paste(element,sep = "")

cluster_list[[element]] <- data_cluster
}

### Making a datatable

matrix_counts <- rbind(data.frame(cluster_list[["cluster1"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster2"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster3"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster4"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster5"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster6"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster7"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster8"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster9"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster10"]])[,c("cluster","counts")])

level_vector <- cluster_vector

### Making violin plot

p <- ggplot(matrix_counts, aes(x=factor(cluster,level_vector), y=log2(counts), color = cluster)) + 
  geom_violin()  + ylim(0,11) + geom_boxplot(width=0.1) + theme_bw() + theme(axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_colour_npg() + geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)

pdf("violinplots-allclusters-H3K27ac-log2.pdf")
print(p)
dev.off()



### Reading data ARID1A

data <- fread("ChIP-seq_HAP1_WT_ARID1A.hg38.raw.counts.txt",header = F)

colnames(data) <- c("chr","pos1","pos2","counts")

data$index <- paste(data$chr,":",data$pos1,"-",data$pos2,sep = "")

### Analysis for ARID1A all regions

### Opening list

cluster_list <- list()

### iterating over clusters

cluster_vector <- c("cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8","cluster9","cluster10")

for (element in cluster_vector){
  
  print(element)
  
  cluster <- read.table(paste("/data/",element,".bed",sep = ""),header = F)
  
  colnames(cluster) <- c("chr","pos1","pos2")
  
  cluster$index <- paste(cluster$chr,":",cluster$pos1,"-",cluster$pos2,sep = "")
  
  data_cluster <- merge(cluster,data,"index")
  
  data_cluster$cluster <- paste(element,sep = "")
  
  cluster_list[[element]] <- data_cluster
}

### Making a datatable

matrix_counts <- rbind(data.frame(cluster_list[["cluster1"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster2"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster3"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster4"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster5"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster6"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster7"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster8"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster9"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster10"]])[,c("cluster","counts")])

level_vector <- cluster_vector

### Making violin plot

p <- ggplot(matrix_counts, aes(x=factor(cluster,level_vector), y=log2(counts), color = cluster)) + 
  geom_violin() + ylim(0,11) + geom_boxplot(width=0.1)  + theme_bw() + theme(axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_colour_npg() + geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)

pdf("violinplots-allclusters-ARID1A-log2.pdf")
print(p)
dev.off()




### Reading data BRD4

data <- fread("ChIP-seq_HAP1_WT_BRD4_r1.hg38.raw.counts.txt",header = F)

colnames(data) <- c("chr","pos1","pos2","counts")

data$index <- paste(data$chr,":",data$pos1,"-",data$pos2,sep = "")

### Analysis for BRD4 all regions

### Opening list

cluster_list <- list()

### Iterating over clusters

cluster_vector <- c("cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8","cluster9","cluster10")

for (element in cluster_vector){
  
  print(element)
  
  cluster <- read.table(paste("/data/",element,".bed",sep = ""),header = F)
  
  colnames(cluster) <- c("chr","pos1","pos2")
  
  cluster$index <- paste(cluster$chr,":",cluster$pos1,"-",cluster$pos2,sep = "")
  
  data_cluster <- merge(cluster,data,"index")
  
  data_cluster$cluster <- paste(element,sep = "")
  
  cluster_list[[element]] <- data_cluster
}

### Making a datatable

matrix_counts <- rbind(data.frame(cluster_list[["cluster1"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster2"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster3"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster4"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster5"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster6"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster7"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster8"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster9"]])[,c("cluster","counts")],data.frame(cluster_list[["cluster10"]])[,c("cluster","counts")])

level_vector <- cluster_vector

### Making violin plots

p <- ggplot(matrix_counts, aes(x=factor(cluster,level_vector), y=log2(counts), color = cluster)) + 
  geom_violin()  + ylim(0,11) + geom_boxplot(width=0.1) + theme_bw() + theme(axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_colour_npg() + geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)

pdf("violinplots-allclusters-BRD4-log2.pdf")
print(p)
dev.off()
