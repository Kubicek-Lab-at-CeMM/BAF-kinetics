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

### atac dtag time-course

FC_atac <- read.table("consensus-FC-atac-data-per-column.txt",header = T)

FC_atac_SMARCA4_tc <- FC_atac[,c("index","log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO",          
                                 "log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO",          
                                 "log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO",          
                                 "log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO",         
                                 "log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO")]

colnames(FC_atac_SMARCA4_tc) <- gsub("log2FoldChange_ATAC.seq_HAP1_","log2FC_dtag_atac_",colnames(FC_atac_SMARCA4_tc))

### chip BRD4 24h and ARID1A

FC_chip <- read.table("macs2-log2-FC-rpm-based.txt",header = T)

FC_chip_BRD4_tc <- FC_chip[,c("index", "FC_WT_SMARCA4dTAG_DMSO_24h_BRD4","FC_WT_SMARCA4dTAG_DMSO_24h_ARID1A")]

### atac ATP depletion time-course + pro seq

FC_atac_ATP <- read.table("ATP-depletion-data-new-deseq2-treatment-per-column-all-regions.txt",header = T)

FC_atac_ATP_tc <- FC_atac_ATP[,c("index_compound","log2FoldChange_N_5min","log2FoldChange_N_10min","log2FoldChange_N_30min","log2FoldChange_N_1h","log2FoldChange_N_6h","log2FoldChange_N_24h")]

merge_file <- read.table("dtag_ATP_overlap_merge.bed",header=F,stringsAsFactors = F)

colnames(merge_file) <- c("chr_dtag","pos1_dtag","pos2_dtag","chr_compound","pos1_compound","pos2_compound")

merge_file$index_dtag <- paste(merge_file$chr_dtag,":",merge_file$pos1_dtag,"-",merge_file$pos2_dtag,sep="")

merge_file$index_compound <- paste(merge_file$chr_compound,":",merge_file$pos1_compound,"-",merge_file$pos2_compound,sep="")

FC_atac_ATP_tc_merge <- merge(FC_atac_ATP_tc,merge_file,"index_compound")

colnames(FC_atac_ATP_tc_merge)[14] <- "index"

FC_pro_seq_ATP <- fread("deseq2-Novartis-proseq-merge-column.txt",data.table = F)

FC_pro_seq_ATP_want <- FC_pro_seq_ATP[,c("index","log2FoldChange_10min_Novartis","log2FoldChange_30min_Novartis","log2FoldChange_1h_Novartis")]

colnames(FC_pro_seq_ATP_want)[1:4] <- c("index_compound",paste(colnames(FC_pro_seq_ATP_want)[2:4],"_proseq",sep=""))

merge_ATP_depl_atac_proseq <- merge(FC_atac_ATP_tc_merge,FC_pro_seq_ATP_want,"index_compound")

merge_ATP_depl_atac_proseq_want <- merge_ATP_depl_atac_proseq[,c("index","log2FoldChange_N_5min","log2FoldChange_N_10min","log2FoldChange_N_30min","log2FoldChange_N_1h","log2FoldChange_N_6h","log2FoldChange_N_24h","log2FoldChange_10min_Novartis_proseq","log2FoldChange_30min_Novartis_proseq","log2FoldChange_1h_Novartis_proseq")]

FC_chip_tc <- fread("log2FC-rpm-based-chip-seq-batch2-dtag-atac-consensus-regions.txt",data.table = F)

### chip SMARCA4 tc

FC_chip_SMARCA4_tc <- FC_chip_tc[,c(1:2,6,10,14)]

### pro seq

FC_pro_seq <- fread("deseq2-3h-dtag-WT-SMARCA4-dtag-vs-3h-DMSO-WT-SMARCA4-dtag.txt",data.table = F)

colnames(FC_pro_seq)[1:6] <- paste(colnames(FC_pro_seq)[1:6],"_proseq_WT_SMARCA4dtag_3h_dtag_vs_DMSO",sep="")

FC_pro_seq_want <- FC_pro_seq[,c("index","log2FoldChange_proseq_WT_SMARCA4dtag_3h_dtag_vs_DMSO")]

### genes

data_gene <- read.csv("BAF_Timecourse_Set3.gene.annotation.singleline.txt",header = T)

colnames(data_gene)[1] <- "index"

data_gene_clean <- data_gene[!data_gene$gene_name == "",]

data_rna <- read.table("rna-seq-all-comp-consensus.txt",header = T)

data_rna_fc <- data_rna[,c(1,3,8,23,38,53,68,83,98,113,128,143,158,173,188,203,218,233, 248, 263, 278, 293, 308, 323, 338, 353, 368, 383, 398, 413, 428, 443, 458, 473, 488, 503, 518, 533, 548, 563, 578, 593)]

write.table(data_rna_fc,"rna-seq-all-comp-consensus-log2FConly.txt",quote = F,sep = "\t",row.names = F)

colnames(data_rna_fc)[2] <- "gene_name"

data_gene_rna <- merge(data_gene_clean,data_rna_fc,"gene_name")

### merging fold change data

data_list_1 <- list(FC_atac_SMARCA4_tc,FC_chip_BRD4_tc,FC_chip_SMARCA4_tc,FC_pro_seq_want)

data_merge_1 <- Reduce(function(x,y) merge(x,y,"index"),data_list_1)

data_merge_2 <- merge_ATP_depl_atac_proseq_want

data_merge_3 <- data_gene_rna[,c("index","log2FoldChange_WT_SMARCA4_F1_dTAG_47_2h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_3h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_6h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_24h_vs_HAP1_WT_SMARCA4_F1_DMSO","log2FoldChange_WT_SMARCA4_F1_dTAG_47_72h_vs_HAP1_WT_SMARCA4_F1_DMSO")]

### looping over clusters

cluster_vector <- c("cluster1","cluster2","cluster3","cluster4","cluster5")

for (element in cluster_vector){
  
  
  cluster = read.table(paste("",element,"-annotated.bed",sep=""),header = F)
  
  colnames(cluster) <- c("index","cluster")
  
  
  
  data_merge_cluster_1 <- merge(data_merge_1,cluster,"index")
  
  data_merge_cluster_2 <- merge(data_merge_2,cluster,"index")
  
  data_merge_cluster_3 <- merge(data_merge_3,cluster,"index")
  
  
  ### calculating median and max/min per sample
  
  median_1 <- apply(data_merge_cluster_1[2:13],2,median_maker)
  
  se_1_pos <- median_1 + apply(data_merge_cluster_1[2:13],2,se_maker)
  
  se_1_neg <- median_1 - (apply(data_merge_cluster_1[2:13],2,se_maker))
  
  
  median_2 <- apply(data_merge_cluster_2[2:10],2,median_maker)
  
  se_2_pos <- median_2 + apply(data_merge_cluster_2[2:10],2,se_maker)
  
  se_2_neg <- median_2 - (apply(data_merge_cluster_2[2:10],2,se_maker))
  
  
  median_3 <- apply(data_merge_cluster_3[2:6],2,median_maker)
  
  se_3_pos <- median_3 + apply(data_merge_cluster_3[2:6],2,se_maker)
  
  se_3_neg <- median_3 - (apply(data_merge_cluster_3[2:6],2,se_maker))
  
  
  median_all <- c(median_1,median_2,median_3)
  
  se_pos_all <- c(se_1_pos,se_2_pos,se_3_pos)
  
  se_neg_all <- c(se_1_neg,se_2_neg,se_3_neg)
  
  
  data <- data.frame(median_all,se_pos_all,se_neg_all)
  
  data$type <- NA
  
  data$timepoint <- NA
  
  data$type[1:5] <- "atac-WT-SMARCA4dtag"
  
  data$type[6] <- "chip-WT-SMARCA4dtag-BRD4"
  
  data$type[7] <- "chip-WT-SMARCA4dtag-ARID1A"
  
  data$type[c(8:11)] <- "chip-WT_SMARCA4dtag-H3K27ac"
  
  data$type[12] <- "proseq-WT_SMARCA4dtag"
  
  data$type[13:18] <- "atac-Novartis-inhibitor"
  
  data$type[19:21] <- "proseq-Novartis-inhibitor"
  
  data$type[22:26] <- "rna-seq-WT-SMARCA4dtag"
  
  data$timepoint[1:26] <- c("2h","3h","6h","24h","72h",replicate(2,"24h"),"3h","6h","24h","72h","3h","5min","10min","30min","1h","6h","24h","10min","30min","1h","2h","3h","6h","24h","72h")
  
  
  # Make the plot
  
  data$timepoint <- factor(data$timepoint, levels=c("5min","10min","30min","1h","2h","3h","6h","24h","72h"))
  
  pdf(paste("lineplot-",element,"-with-all-data.pdf",sep=""))
  print(ggplot(data=data, aes(x=timepoint, y=median_all, group=type)) +
          geom_line(aes(color=type)) +
          geom_point(aes(color=type), size=2) +
          geom_ribbon(aes(ymin = se_neg_all, ymax = se_pos_all, fill = type),alpha = 0.3) +
          ylab("Fold change") + theme_bw() + ylim(-2.5,2.5) + theme(axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
  dev.off()
  
}








