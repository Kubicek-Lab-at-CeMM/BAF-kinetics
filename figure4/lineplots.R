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

FC_atac_interest <- FC_atac[,c("index",
                               "log2FoldChange_ATAC.seq_HAP1_A2_A4_F5_dTAG_47_2h_vs_DMSO",
                               "log2FoldChange_ATAC.seq_HAP1_A2_A4_F5_dTAG_47_3h_vs_DMSO",               
                               "log2FoldChange_ATAC.seq_HAP1_A2_A4_F5_dTAG_47_6h_vs_DMSO",
                               "log2FoldChange_ATAC.seq_HAP1_A2_A4_F5_dTAG_47_24h_vs_DMSO",              
                               "log2FoldChange_ATAC.seq_HAP1_A2_A4_F5_dTAG_47_72h_vs_DMSO",
                               "log2FoldChange_ATAC.seq_HAP1_BRG_BRM_C8_dTAG_47_2h_vs_DMSO",
                               "log2FoldChange_ATAC.seq_HAP1_BRG_BRM_C8_dTAG_47_3h_vs_DMSO",             
                               "log2FoldChange_ATAC.seq_HAP1_BRG_BRM_C8_dTAG_47_6h_vs_DMSO",
                               "log2FoldChange_ATAC.seq_HAP1_BRG_BRM_C8_dTAG_47_24h_vs_DMSO",            
                               "log2FoldChange_ATAC.seq_HAP1_BRG_BRM_C8_dTAG_47_72h_vs_DMSO",
                               "log2FoldChange_ATAC.seq_HAP1_CC1_CC2_F1_dTAG_47_3h_vs_combined_control", 
                               "log2FoldChange_ATAC.seq_HAP1_CC1_CC2_F1_dTAG_47_6h_vs_combined_control",
                               "log2FoldChange_ATAC.seq_HAP1_CC1_CC2_F1_dTAG_47_24h_vs_combined_control",
                               "log2FoldChange_ATAC.seq_HAP1_CC1_CC2_F1_dTAG_47_72h_vs_combined_control")]

colnames(FC_atac_interest) <- gsub("log2FoldChange_ATAC.seq_HAP1_","log2FC_dtag_atac_",colnames(FC_atac_interest))

### chip BRD4 24h and ARID1A and H3K27ac

FC_chip <- read.table("macs2-log2-FC-rpm-based.txt",header = T)

FC_chip_tc_old <- FC_chip[,c("index", "FC_WT_SMARCA4dTAG_DMSO_24h_BRD4","FC_WT_SMARCA4dTAG_DMSO_24h_ARID1A")]

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

### chip batch2

FC_chip_SL_tc <- FC_chip_tc[,c(1,18,22,26,30)]

### genes

data_gene <- read.csv("BAF_Timecourse_Set3.gene.annotation.singleline.txt",header = T)

colnames(data_gene)[1] <- "index"

data_gene_clean <- data_gene[!data_gene$gene_name == "",]

data_rna <- read.table("rna-seq-all-comp-consensus-diff-timecourse.txt",header = T)

data_rna_fc <- data_rna[,c(1,3,8,23,38,53,68,83,98,113,128,143,158,173,188,203,218,233, 248, 263, 278, 293, 308, 323, 338, 353, 368, 383, 398, 413, 428, 443, 458, 473, 488, 503, 518, 533, 548, 563, 578, 593)]

colnames(data_rna_fc)[2] <- "gene_name"

data_gene_rna <- merge(data_gene_clean,data_rna_fc,"gene_name")


### merging fold change data

data_list_1 <- list(FC_atac_interest,FC_chip_tc_old,FC_chip_SL_tc)

data_merge_1 <- Reduce(function(x,y) merge(x,y,"index"),data_list_1)

data_merge_2 <- merge_ATP_depl_atac_proseq_want

data_merge_3 <- data_gene_rna[,c("index","gene_name","log2FoldChange_A2_A4_F5_dTAG_47_24h_vs_HAP1_A2_A4_F5_DMSO",                                            
                                "log2FoldChange_A2_A4_F5_dTAG_47_2h_vs_HAP1_A2_A4_F5_DMSO",                                               
                                "log2FoldChange_A2_A4_F5_dTAG_47_3h_vs_HAP1_A2_A4_F5_DMSO",                                              
                                "log2FoldChange_A2_A4_F5_dTAG_47_6h_vs_HAP1_A2_A4_F5_DMSO",                                               
                                "log2FoldChange_A2_A4_F5_dTAG_47_72h_vs_HAP1_A2_A4_F5_DMSO",
                                "log2FoldChange_CC1_CC2_F1_dTAG_47_24h_vs_HAP1_CC1_CC2_F1_DMSO",                                          
                                "log2FoldChange_CC1_CC2_F1_dTAG_47_2h_vs_HAP1_CC1_CC2_F1_DMSO",                                           
                                "log2FoldChange_CC1_CC2_F1_dTAG_47_3h_vs_HAP1_CC1_CC2_F1_DMSO",                                           
                                "log2FoldChange_CC1_CC2_F1_dTAG_47_6h_vs_HAP1_CC1_CC2_F1_DMSO",                                           
                                "log2FoldChange_CC1_CC2_F1_dTAG_47_72h_vs_HAP1_CC1_CC2_F1_DMSO",
                                "log2FoldChange_BRG_BRM_C8_dTAG_47_24h_vs_HAP1_BRG_BRM_C8_DMSO",                                          
                                "log2FoldChange_BRG_BRM_C8_dTAG_47_2h_vs_HAP1_BRG_BRM_C8_DMSO",                                           
                                "log2FoldChange_BRG_BRM_C8_dTAG_47_3h_vs_HAP1_BRG_BRM_C8_DMSO",                                           
                                "log2FoldChange_BRG_BRM_C8_dTAG_47_6h_vs_HAP1_BRG_BRM_C8_DMSO",                                           
                                "log2FoldChange_BRG_BRM_C8_dTAG_47_72h_vs_HAP1_BRG_BRM_C8_DMSO"
                                )]


### looping over clusters

cluster_vector <- c("cluster6","cluster7","cluster8","cluster9","cluster10","cluster11")

for (element in cluster_vector){
  
  
  cluster = read.table(paste("",element,"-annotated.bed",sep=""),header = F)
  
  colnames(cluster) <- c("index","cluster")
  
  
  
  data_merge_cluster_1 <- merge(data_merge_1,cluster,"index")
  
  data_merge_cluster_2 <- merge(data_merge_2,cluster,"index")
  
  print(length(data_merge_cluster_2$index))
  
  data_merge_cluster_3 <- merge(data_merge_3,cluster,"index")
  
  data_gene_ARID2KO <- read.table("ARID2KO_SMARCA4dtag_genes_differential_unique.txt",header = F)
  
  colnames(data_gene_ARID2KO) <- "gene_name"
  
  data_merge_cluster3_ARID2KO <- merge(data_gene_ARID2KO,data_merge_cluster_3,"gene_name")
  
  data_gene_BRG <- read.table("BRG_BRMdtag_genes_differential_unique.txt",header = F)
  
  colnames(data_gene_BRG) <- "gene_name"
  
  data_merge_cluster3_BRG <- merge(data_gene_BRG,data_merge_cluster_3,"gene_name")
  
  data_gene_CC1 <- read.table("CC1_CC2dtag_genes_differential_unique.txt",header = F)
  
  colnames(data_gene_CC1) <- "gene_name"
  
  data_merge_cluster3_CC1 <- merge(data_gene_CC1,data_merge_cluster_3,"gene_name")
  
  
  ### calculating median and max/min per sample
  
  median_1 <- apply(data_merge_cluster_1[2:21],2,median_maker)
  
  se_1_pos <- median_1 + apply(data_merge_cluster_1[2:21],2,se_maker)
  
  se_1_neg <- median_1 - (apply(data_merge_cluster_1[2:21],2,se_maker))
  
  
  median_2 <- apply(data_merge_cluster_2[2:10],2,median_maker)
  
  se_2_pos <- median_2 + apply(data_merge_cluster_2[2:10],2,se_maker)
  
  se_2_neg <- median_2 - (apply(data_merge_cluster_2[2:10],2,se_maker))
  
  
  median_3 <- apply(data_merge_cluster3_ARID2KO[3:7],2,median_maker)
  
  se_3_pos <- median_3 + apply(data_merge_cluster3_ARID2KO[3:7],2,se_maker)
  
  se_3_neg <- median_3 - (apply(data_merge_cluster3_ARID2KO[3:7],2,se_maker))
  
  
  median_4 <- apply(data_merge_cluster3_BRG[13:17],2,median_maker)
  
  se_4_pos <- median_4 + apply(data_merge_cluster3_BRG[13:17],2,se_maker)
  
  se_4_neg <- median_4 - (apply(data_merge_cluster3_BRG[13:17],2,se_maker))
  
  
  median_5 <- apply(data_merge_cluster3_CC1[8:12],2,median_maker)
  
  se_5_pos <- median_5 + apply(data_merge_cluster3_CC1[8:12],2,se_maker)
  
  se_5_neg <- median_5 - (apply(data_merge_cluster3_CC1[8:12],2,se_maker))
  
  
  median_all <- c(median_1,median_2,median_3,median_4,median_5)
  
  se_pos_all <- c(se_1_pos,se_2_pos,se_3_pos,se_4_pos,se_5_pos)
  
  se_neg_all <- c(se_1_neg,se_2_neg,se_3_neg,se_4_neg,se_5_neg)
  
  
  data <- data.frame(median_all,se_pos_all,se_neg_all)
  
  data$type <- NA
  
  data$timepoint <- NA
  
  data$type[1:5] <- "atac-ARID2KO-SMARCA4dtag"
  
  data$type[6:10] <- "atac-BRGKO-BRMdtag"
  
  data$type[11:14] <- "atac-CC1KO-CC2dtag"
  
  data$type[15] <- "chip-WT-SMARCA4dtag-BRD4"
  
  data$type[16] <- "chip-WT-SMARCA4dtag-ARID1A"
  
  data$type[17:20] <- "chip-SMARCA4KO-SMARCA2dtag-H3K27ac"
  
  data$type[21:26] <- "atac-Novartis-inhibitor"
  
  data$type[27:29] <- "proseq-Novartis-inhibitor"
  
  data$type[30:34] <- "rna-seq-ARID2KO-SMARCA4dtag"
  
  data$type[35:39] <- "rna-seq-BRG-BRMdtag"
  
  data$type[40:44] <- "rna-seq-CC1-CC2dtag"
  
  data$timepoint[1:44] <- c("2h","3h","6h","24h","72h","2h","3h","6h","24h","72h","3h","6h","24h","72h",replicate(2,"24h"),"3h","6h","24h","72h","5min","10min","30min","1h","6h","24h","10min","30min","1h","24h","2h","3h","6h","72h","24h","2h","3h","6h","72h","24h","2h","3h","6h","72h")
  
  
  # Make the plot
  
  data$timepoint <- factor(data$timepoint, levels=c("5min","10min","30min","1h","2h","3h","6h","24h","72h"))
  
  pdf(paste("lineplot-",element,"-with-all-data-diff-rnaseq-per-cluster.pdf",sep=""))
  print(ggplot(data=data, aes(x=timepoint, y=median_all, group=type)) +
          geom_line(aes(color=type)) +
          geom_point(aes(color=type), size=2) +
          geom_ribbon(aes(ymin = se_neg_all, ymax = se_pos_all, fill = type),alpha = 0.3) +
          ylab("Fold change") + theme_bw() + ylim(-2.5,2.5) + theme(axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
  dev.off()
  
}
