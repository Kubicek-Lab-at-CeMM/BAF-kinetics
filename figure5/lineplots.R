### getting packages

library(ggplot2)
library(reshape2)
library(ggsci)

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



#### super enhancer

### reading in FC data

data <- read.table("supenhancer-diff-any-condition-merged-deseq2-output.txt",header = T)

sup_enh <- data.frame(unique(data$index))

colnames(sup_enh) <- "index"

data_chip_atac <- read.table("log2FC-rpm-based-chip-seq-batch2-dtag-atac-consensus-regions.txt",header = T)

data_merge <- merge(sup_enh,data_chip_atac,"index")

FC_atac <- read.table("consensus-FC-atac-data-per-column.txt",header = T)

FC_atac_interest <- FC_atac[,c("index",
                               "log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO",          
                               "log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO",          
                               "log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO",          
                               "log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO",         
                               "log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO",         
                               "log2FoldChange_ATAC.seq_HAP1_BRGKO_vs_WT",
                               "log2FoldChange_ATAC.seq_HAP1_BRG_BRM_C8_dTAG_47_2h_vs_DMSO",
                               "log2FoldChange_ATAC.seq_HAP1_BRG_BRM_C8_dTAG_47_3h_vs_DMSO",             
                               "log2FoldChange_ATAC.seq_HAP1_BRG_BRM_C8_dTAG_47_6h_vs_DMSO",
                               "log2FoldChange_ATAC.seq_HAP1_BRG_BRM_C8_dTAG_47_24h_vs_DMSO",            
                               "log2FoldChange_ATAC.seq_HAP1_BRG_BRM_C8_dTAG_47_72h_vs_DMSO",
                               "log2FoldChange_ATAC.seq_HAP1_BRMKO_vs_WT",
                               "log2FoldChange_ATAC.seq_HAP1_ARID2KO_vs_WT",                             
                               "log2FoldChange_ATAC.seq_HAP1_A2_A4_F5_dTAG_47_2h_vs_DMSO",               
                               "log2FoldChange_ATAC.seq_HAP1_A2_A4_F5_dTAG_47_3h_vs_DMSO",               
                               "log2FoldChange_ATAC.seq_HAP1_A2_A4_F5_dTAG_47_6h_vs_DMSO",               
                               "log2FoldChange_ATAC.seq_HAP1_A2_A4_F5_dTAG_47_24h_vs_DMSO",              
                               "log2FoldChange_ATAC.seq_HAP1_A2_A4_F5_dTAG_47_72h_vs_DMSO",
                               "log2FoldChange_ATAC.seq_HAP1_CC1KO_vs_WT",
                               "log2FoldChange_ATAC.seq_HAP1_CC2_vs_WT",
                               "log2FoldChange_ATAC.seq_HAP1_CC1_CC2_F1_dTAG_47_3h_vs_combined_control", 
                               "log2FoldChange_ATAC.seq_HAP1_CC1_CC2_F1_dTAG_47_6h_vs_combined_control", 
                               "log2FoldChange_ATAC.seq_HAP1_CC1_CC2_F1_dTAG_47_24h_vs_combined_control",
                               "log2FoldChange_ATAC.seq_HAP1_CC1_CC2_F1_dTAG_47_72h_vs_combined_control"
                               
)]

colnames(FC_atac_interest) <- gsub("log2FoldChange_ATAC.seq_HAP1_","log2FC_dtag_atac_",colnames(FC_atac_interest))

data_merge_all <- merge(data_merge,FC_atac_interest,"index")

### calculating median per sample

median_all <- data.frame(apply(data_merge_all[2:length(data_merge_all)],2,median_maker))

colnames(median_all) <- "median"

se_all_pos <- data.frame(median_all$median + apply(data_merge_all[2:length(data_merge_all)],2,se_maker))

colnames(se_all_pos) <- "se_pos"

se_all_neg <- data.frame(median_all$median - apply(data_merge_all[2:length(data_merge_all)],2,se_maker))

colnames(se_all_neg) <- "se_neg"

### lineplots with KO one with standard error

median_all_interest <- data.frame(median_all[c("log2FCWT_SMARCA4dTAG_3h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_6h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_24h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_72h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA2KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_3h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_6h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_24h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_72h_H3K27ac_vs_DMSO",
                                               "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO",
                                               "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO",  
                                               "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO",  
                                               "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO",
                                               "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO", 
                                               "log2FC_dtag_atac_BRGKO_vs_WT",
                                               "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_2h_vs_DMSO",
                                               "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_3h_vs_DMSO",   
                                               "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_6h_vs_DMSO",    
                                               "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_24h_vs_DMSO",  
                                               "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_72h_vs_DMSO", 
                                               "log2FC_dtag_atac_BRMKO_vs_WT",
                                               "log2FC_dtag_atac_ARID2KO_vs_WT",
                                               "log2FC_dtag_atac_A2_A4_F5_dTAG_47_2h_vs_DMSO",
                                               "log2FC_dtag_atac_A2_A4_F5_dTAG_47_3h_vs_DMSO",
                                               "log2FC_dtag_atac_A2_A4_F5_dTAG_47_6h_vs_DMSO",
                                               "log2FC_dtag_atac_A2_A4_F5_dTAG_47_24h_vs_DMSO",
                                               "log2FC_dtag_atac_A2_A4_F5_dTAG_47_72h_vs_DMSO",
                                               "log2FC_dtag_atac_CC1KO_vs_WT",
                                               "log2FC_dtag_atac_CC2_vs_WT",
                                               "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_3h_vs_combined_control",
                                               "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_6h_vs_combined_control",
                                               "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_24h_vs_combined_control",
                                               "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_72h_vs_combined_control"
                                               ),])

colnames(median_all_interest) <- "median"

median_all_interest$sample <- c("log2FCWT_SMARCA4dTAG_3h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_6h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_24h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_72h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA2KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_3h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_6h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_24h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_72h_H3K27ac_vs_DMSO",
                                "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO",
                                "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO",  
                                "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO",  
                                "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO",
                                "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO", 
                                "log2FC_dtag_atac_BRGKO_vs_WT",
                                "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_2h_vs_DMSO",
                                "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_3h_vs_DMSO",   
                                "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_6h_vs_DMSO",    
                                "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_24h_vs_DMSO",  
                                "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_72h_vs_DMSO", 
                                "log2FC_dtag_atac_BRMKO_vs_WT",
                                "log2FC_dtag_atac_ARID2KO_vs_WT",
                                "log2FC_dtag_atac_A2_A4_F5_dTAG_47_2h_vs_DMSO",
                                "log2FC_dtag_atac_A2_A4_F5_dTAG_47_3h_vs_DMSO",
                                "log2FC_dtag_atac_A2_A4_F5_dTAG_47_6h_vs_DMSO",
                                "log2FC_dtag_atac_A2_A4_F5_dTAG_47_24h_vs_DMSO",
                                "log2FC_dtag_atac_A2_A4_F5_dTAG_47_72h_vs_DMSO",
                                "log2FC_dtag_atac_CC1KO_vs_WT",
                                "log2FC_dtag_atac_CC2_vs_WT",
                                "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_3h_vs_combined_control",
                                "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_6h_vs_combined_control",
                                "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_24h_vs_combined_control",
                                "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_72h_vs_combined_control")

median_all_interest$timepoint <- c("3h","6h","24h","72h","SMARCA4KO","SMARCA2KO","3h","6h","24h","72h","2h","3h","6h","24h","72h","SMARCA4KO","2h","3h","6h","24h","72h","SMARCA2KO","ARID2KO","2h","3h","6h","24h","72h","CC1KO","CC2KO","3h","6h","24h","72h")

median_all_interest$type <- c(replicate(5,"WT_SMARCA4dtag_H3K27ac_chip"),replicate(5,"BRG_BRMdtag_H3K27ac_chip"),replicate(6,"WT_SMARCA4dtag_atac"),replicate(6,"BRG_BRMdtag_atac"),replicate(6,"ARID2KO_SMARCA4dtag_atac"),replicate(6,"CC1KO_CC2dtag_atac"))

median_all_interest$se_pos <- se_all_pos[c("log2FCWT_SMARCA4dTAG_3h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_6h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_24h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_72h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA2KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_3h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_6h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_24h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_72h_H3K27ac_vs_DMSO",
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO",
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO",  
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO",  
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO",
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO", 
                                           "log2FC_dtag_atac_BRGKO_vs_WT",
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_2h_vs_DMSO",
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_3h_vs_DMSO",   
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_6h_vs_DMSO",    
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_24h_vs_DMSO",  
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_72h_vs_DMSO", 
                                           "log2FC_dtag_atac_BRMKO_vs_WT",
                                           "log2FC_dtag_atac_ARID2KO_vs_WT",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_2h_vs_DMSO",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_3h_vs_DMSO",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_6h_vs_DMSO",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_24h_vs_DMSO",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_72h_vs_DMSO",
                                           "log2FC_dtag_atac_CC1KO_vs_WT",
                                           "log2FC_dtag_atac_CC2_vs_WT",
                                           "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_3h_vs_combined_control",
                                           "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_6h_vs_combined_control",
                                           "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_24h_vs_combined_control",
                                           "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_72h_vs_combined_control"),]


median_all_interest$se_neg <- se_all_neg[c("log2FCWT_SMARCA4dTAG_3h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_6h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_24h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_72h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA2KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_3h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_6h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_24h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_72h_H3K27ac_vs_DMSO",
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO",
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO",  
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO",  
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO",
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO", 
                                           "log2FC_dtag_atac_BRGKO_vs_WT",
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_2h_vs_DMSO",
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_3h_vs_DMSO",   
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_6h_vs_DMSO",    
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_24h_vs_DMSO",  
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_72h_vs_DMSO", 
                                           "log2FC_dtag_atac_BRMKO_vs_WT",
                                           "log2FC_dtag_atac_ARID2KO_vs_WT",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_2h_vs_DMSO",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_3h_vs_DMSO",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_6h_vs_DMSO",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_24h_vs_DMSO",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_72h_vs_DMSO",
                                           "log2FC_dtag_atac_CC1KO_vs_WT",
                                           "log2FC_dtag_atac_CC2_vs_WT",
                                           "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_3h_vs_combined_control",
                                           "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_6h_vs_combined_control",
                                           "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_24h_vs_combined_control",
                                           "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_72h_vs_combined_control"),]


median_all_interest$annotation <- "superenhancer"





#### enhancer

### reading in FC data

data <- read.table("enhancer-diff-wo-supenhancer-all.bed",header = F)

data$index <- paste(data$V1,":",data$V2,"-",data$V3,sep = "")

enh <- data.frame(unique(data$index))

colnames(enh) <- "index"

data_chip_atac <- read.table("log2FC-rpm-based-chip-seq-batch2-dtag-atac-consensus-regions.txt",header = T)

data_merge <- merge(enh,data_chip_atac,"index")

FC_atac <- read.table("consensus-FC-atac-data-per-column.txt",header = T)

FC_atac_interest <- FC_atac[,c("index",
                               "log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO",          
                               "log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO",          
                               "log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO",          
                               "log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO",         
                               "log2FoldChange_ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO",         
                               "log2FoldChange_ATAC.seq_HAP1_BRGKO_vs_WT",
                               "log2FoldChange_ATAC.seq_HAP1_BRG_BRM_C8_dTAG_47_2h_vs_DMSO",
                               "log2FoldChange_ATAC.seq_HAP1_BRG_BRM_C8_dTAG_47_3h_vs_DMSO",             
                               "log2FoldChange_ATAC.seq_HAP1_BRG_BRM_C8_dTAG_47_6h_vs_DMSO",
                               "log2FoldChange_ATAC.seq_HAP1_BRG_BRM_C8_dTAG_47_24h_vs_DMSO",            
                               "log2FoldChange_ATAC.seq_HAP1_BRG_BRM_C8_dTAG_47_72h_vs_DMSO",
                               "log2FoldChange_ATAC.seq_HAP1_BRMKO_vs_WT",
                               "log2FoldChange_ATAC.seq_HAP1_ARID2KO_vs_WT",                             
                               "log2FoldChange_ATAC.seq_HAP1_A2_A4_F5_dTAG_47_2h_vs_DMSO",               
                               "log2FoldChange_ATAC.seq_HAP1_A2_A4_F5_dTAG_47_3h_vs_DMSO",               
                               "log2FoldChange_ATAC.seq_HAP1_A2_A4_F5_dTAG_47_6h_vs_DMSO",               
                               "log2FoldChange_ATAC.seq_HAP1_A2_A4_F5_dTAG_47_24h_vs_DMSO",              
                               "log2FoldChange_ATAC.seq_HAP1_A2_A4_F5_dTAG_47_72h_vs_DMSO",
                               "log2FoldChange_ATAC.seq_HAP1_CC1KO_vs_WT",
                               "log2FoldChange_ATAC.seq_HAP1_CC2_vs_WT",
                               "log2FoldChange_ATAC.seq_HAP1_CC1_CC2_F1_dTAG_47_3h_vs_combined_control", 
                               "log2FoldChange_ATAC.seq_HAP1_CC1_CC2_F1_dTAG_47_6h_vs_combined_control", 
                               "log2FoldChange_ATAC.seq_HAP1_CC1_CC2_F1_dTAG_47_24h_vs_combined_control",
                               "log2FoldChange_ATAC.seq_HAP1_CC1_CC2_F1_dTAG_47_72h_vs_combined_control"
                               
)]

colnames(FC_atac_interest) <- gsub("log2FoldChange_ATAC.seq_HAP1_","log2FC_dtag_atac_",colnames(FC_atac_interest))

data_merge_all <- merge(data_merge,FC_atac_interest,"index")

### calculating median per sample

median_all <- data.frame(apply(data_merge_all[2:length(data_merge_all)],2,median_maker))

colnames(median_all) <- "median"

se_all_pos <- data.frame(median_all$median + apply(data_merge_all[2:length(data_merge_all)],2,se_maker))

colnames(se_all_pos) <- "se_pos"

se_all_neg <- data.frame(median_all$median - apply(data_merge_all[2:length(data_merge_all)],2,se_maker))

colnames(se_all_neg) <- "se_neg"

### lineplots with KO one with standard error

median_all_interest_enh <- data.frame(median_all[c("log2FCWT_SMARCA4dTAG_3h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_6h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_24h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_72h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA2KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_3h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_6h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_24h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_72h_H3K27ac_vs_DMSO",
                                               "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO",
                                               "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO",  
                                               "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO",  
                                               "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO",
                                               "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO", 
                                               "log2FC_dtag_atac_BRGKO_vs_WT",
                                               "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_2h_vs_DMSO",
                                               "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_3h_vs_DMSO",   
                                               "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_6h_vs_DMSO",    
                                               "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_24h_vs_DMSO",  
                                               "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_72h_vs_DMSO", 
                                               "log2FC_dtag_atac_BRMKO_vs_WT",
                                               "log2FC_dtag_atac_ARID2KO_vs_WT",
                                               "log2FC_dtag_atac_A2_A4_F5_dTAG_47_2h_vs_DMSO",
                                               "log2FC_dtag_atac_A2_A4_F5_dTAG_47_3h_vs_DMSO",
                                               "log2FC_dtag_atac_A2_A4_F5_dTAG_47_6h_vs_DMSO",
                                               "log2FC_dtag_atac_A2_A4_F5_dTAG_47_24h_vs_DMSO",
                                               "log2FC_dtag_atac_A2_A4_F5_dTAG_47_72h_vs_DMSO",
                                               "log2FC_dtag_atac_CC1KO_vs_WT",
                                               "log2FC_dtag_atac_CC2_vs_WT",
                                               "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_3h_vs_combined_control",
                                               "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_6h_vs_combined_control",
                                               "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_24h_vs_combined_control",
                                               "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_72h_vs_combined_control"
),])

colnames(median_all_interest_enh) <- "median"

median_all_interest_enh$sample <- c("log2FCWT_SMARCA4dTAG_3h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_6h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_24h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_72h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA2KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_3h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_6h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_24h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_72h_H3K27ac_vs_DMSO",
                                "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO",
                                "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO",  
                                "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO",  
                                "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO",
                                "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO", 
                                "log2FC_dtag_atac_BRGKO_vs_WT",
                                "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_2h_vs_DMSO",
                                "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_3h_vs_DMSO",   
                                "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_6h_vs_DMSO",    
                                "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_24h_vs_DMSO",  
                                "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_72h_vs_DMSO", 
                                "log2FC_dtag_atac_BRMKO_vs_WT",
                                "log2FC_dtag_atac_ARID2KO_vs_WT",
                                "log2FC_dtag_atac_A2_A4_F5_dTAG_47_2h_vs_DMSO",
                                "log2FC_dtag_atac_A2_A4_F5_dTAG_47_3h_vs_DMSO",
                                "log2FC_dtag_atac_A2_A4_F5_dTAG_47_6h_vs_DMSO",
                                "log2FC_dtag_atac_A2_A4_F5_dTAG_47_24h_vs_DMSO",
                                "log2FC_dtag_atac_A2_A4_F5_dTAG_47_72h_vs_DMSO",
                                "log2FC_dtag_atac_CC1KO_vs_WT",
                                "log2FC_dtag_atac_CC2_vs_WT",
                                "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_3h_vs_combined_control",
                                "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_6h_vs_combined_control",
                                "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_24h_vs_combined_control",
                                "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_72h_vs_combined_control")

median_all_interest_enh$timepoint <- c("3h","6h","24h","72h","SMARCA4KO","SMARCA2KO","3h","6h","24h","72h","2h","3h","6h","24h","72h","SMARCA4KO","2h","3h","6h","24h","72h","SMARCA2KO","ARID2KO","2h","3h","6h","24h","72h","CC1KO","CC2KO","3h","6h","24h","72h")

median_all_interest_enh$type <- c(replicate(5,"WT_SMARCA4dtag_H3K27ac_chip"),replicate(5,"BRG_BRMdtag_H3K27ac_chip"),replicate(6,"WT_SMARCA4dtag_atac"),replicate(6,"BRG_BRMdtag_atac"),replicate(6,"ARID2KO_SMARCA4dtag_atac"),replicate(6,"CC1KO_CC2dtag_atac"))

median_all_interest_enh$se_pos <- se_all_pos[c("log2FCWT_SMARCA4dTAG_3h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_6h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_24h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_72h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA2KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_3h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_6h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_24h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_72h_H3K27ac_vs_DMSO",
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO",
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO",  
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO",  
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO",
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO", 
                                           "log2FC_dtag_atac_BRGKO_vs_WT",
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_2h_vs_DMSO",
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_3h_vs_DMSO",   
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_6h_vs_DMSO",    
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_24h_vs_DMSO",  
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_72h_vs_DMSO", 
                                           "log2FC_dtag_atac_BRMKO_vs_WT",
                                           "log2FC_dtag_atac_ARID2KO_vs_WT",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_2h_vs_DMSO",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_3h_vs_DMSO",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_6h_vs_DMSO",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_24h_vs_DMSO",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_72h_vs_DMSO",
                                           "log2FC_dtag_atac_CC1KO_vs_WT",
                                           "log2FC_dtag_atac_CC2_vs_WT",
                                           "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_3h_vs_combined_control",
                                           "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_6h_vs_combined_control",
                                           "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_24h_vs_combined_control",
                                           "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_72h_vs_combined_control"),]


median_all_interest_enh$se_neg <- se_all_neg[c("log2FCWT_SMARCA4dTAG_3h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_6h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_24h_H3K27ac_vs_DMSO","log2FCWT_SMARCA4dTAG_72h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA2KO_DMSO_H3K27ac_vs_WT_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_3h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_6h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_24h_H3K27ac_vs_DMSO","log2FCSMARCA4KO_SMARCA2dTAG_72h_H3K27ac_vs_DMSO",
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO",
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO",  
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO",  
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO",
                                           "log2FC_dtag_atac_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO", 
                                           "log2FC_dtag_atac_BRGKO_vs_WT",
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_2h_vs_DMSO",
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_3h_vs_DMSO",   
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_6h_vs_DMSO",    
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_24h_vs_DMSO",  
                                           "log2FC_dtag_atac_BRG_BRM_C8_dTAG_47_72h_vs_DMSO", 
                                           "log2FC_dtag_atac_BRMKO_vs_WT",
                                           "log2FC_dtag_atac_ARID2KO_vs_WT",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_2h_vs_DMSO",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_3h_vs_DMSO",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_6h_vs_DMSO",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_24h_vs_DMSO",
                                           "log2FC_dtag_atac_A2_A4_F5_dTAG_47_72h_vs_DMSO",
                                           "log2FC_dtag_atac_CC1KO_vs_WT",
                                           "log2FC_dtag_atac_CC2_vs_WT",
                                           "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_3h_vs_combined_control",
                                           "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_6h_vs_combined_control",
                                           "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_24h_vs_combined_control",
                                           "log2FC_dtag_atac_CC1_CC2_F1_dTAG_47_72h_vs_combined_control"),]



median_all_interest_enh$annotation <- "enhancer"



### merging lists

data_plot <- rbind(median_all_interest,median_all_interest_enh)

data_plot$line <- paste(data_plot$type,"_",data_plot$annotation,sep = "")

data_plot[c(69:72),] <- data_plot[c(5,16,39,50),]

data_plot$line[69:72] <- c("BRG_BRMdtag_H3K27ac_chip_superenhancer","BRG_BRMdtag_atac_superenhancer","BRG_BRMdtag_H3K27ac_chip_enhancer","BRG_BRMdtag_atac_enhancer")

data_plot_WT <- data_plot[data_plot$line %in% c("WT_SMARCA4dtag_H3K27ac_chip_superenhancer","WT_SMARCA4dtag_atac_superenhancer","WT_SMARCA4dtag_atac_enhancer","WT_SMARCA4dtag_H3K27ac_chip_enhancer"),]

data_plot_WT$timepoint <- factor(data_plot_WT$timepoint, levels=c("2h","3h","6h","24h","72h","SMARCA4KO"))


data_plot_BRGBRM <- data_plot[data_plot$line %in% c("BRG_BRMdtag_H3K27ac_chip_superenhancer","BRG_BRMdtag_atac_superenhancer","BRG_BRMdtag_H3K27ac_chip_enhancer","BRG_BRMdtag_atac_enhancer"),]

data_plot_BRGBRM$timepoint <- factor(data_plot_BRGBRM$timepoint, levels=c("SMARCA4KO","SMARCA2KO","2h","3h","6h","24h","72h"))



data_plot_ARID2KOSMARCA4dtag <- data_plot[data_plot$line %in% c("ARID2KO_SMARCA4dtag_atac_superenhancer","ARID2KO_SMARCA4dtag_atac_enhancer"),]

data_plot_ARID2KOSMARCA4dtag$timepoint <- factor(data_plot_ARID2KOSMARCA4dtag$timepoint, levels=c("ARID2KO","2h","3h","6h","24h","72h"))


data_plot_CC1CC2dtag <- data_plot[data_plot$line %in% c("CC1KO_CC2dtag_atac_superenhancer","CC1KO_CC2dtag_atac_enhancer"),]

data_plot_CC1CC2dtag$timepoint <- factor(data_plot_CC1CC2dtag$timepoint, levels=c("CC1KO","CC2KO","3h","6h","24h","72h"))


pdf(paste("lineplot-diff-supenhancer-regions-enhancer-regions-log2FC-atac-timecourse-chip-H3K27ac-WT-SMARCA4dtag.pdf",sep=""))
print(ggplot(data=data_plot_WT, aes(x=timepoint, y=median, group=line)) +
        geom_line(aes(color=line)) +
        geom_point(aes(color=line), size=2) +
        geom_ribbon(aes(ymin = se_neg, ymax = se_pos, fill = line),alpha = 0.3) + ylim(-1.4,0.7) +
        ylab("median log2 FC") + theme_bw() + theme(axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
dev.off()




data_plot_BRGBRM$timepoint <- factor(data_plot_BRGBRM$timepoint, levels=c("SMARCA2KO","2h","3h","6h","24h","72h","SMARCA4KO"))

pdf(paste("lineplot-diff-supenhancer-regions-enhancer-regions-log2FC-atac-timecourse-chip-H3K27ac-BRG-BRMdtag-wSMARCA4end.pdf",sep=""))
print(ggplot(data=data_plot_BRGBRM, aes(x=timepoint, y=median, group=line)) +
        geom_line(aes(color=line)) +
        geom_point(aes(color=line), size=2) +
        geom_ribbon(aes(ymin = se_neg, ymax = se_pos, fill = line),alpha = 0.3) + ylim(-1.4,0.7) +
        ylab("median log2 FC") + theme_bw() + theme(axis.text.x = element_text(angle = 90),panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
dev.off()

