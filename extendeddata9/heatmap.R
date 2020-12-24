### getting packages

library(data.table)
library(factoextra)
library(NbClust)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(plyr)
library(ggplot2)
library(hexbin)
library(viridis)
library(scales)
library(reshape)
library(ggsci)

### setting working directory

setwd("../data/")

### reading macs2 qn count data chip batch 1

data_macs2_qn_batch1 <- read.csv("Chip.matrix_qn.norm.for.atac.consensus.batch3.regions.csv")

colnames(data_macs2_qn_batch1)[1] <- "index"

columns_wanted <- c("index","WT_SMARCA4dTAG_DMSO_ARID1A",
                    "WT_SMARCA4dTAG_dTAG47_24h_ARID1A",
                    "WT_SMARCA4dTAG_DMSO_BRD4",
                    "WT_SMARCA4dTAG_dTAG47_24h_BRD4",
                    "WT_SMARCA4dTAG_DMSO_H3K27ac",
                    "WT_SMARCA4dTAG_dTAG47_24h_H3K27ac",
                    "SMARCA4KO_SMARCA2dTAG_DMSO_ARID1A",
                    "SMARCA4KO_SMARCA2dTAG_dTAG47_24h_ARID1A",
                    "SMARCA4KO_SMARCA2dTAG_DMSO_BRD4",
                    "SMARCA4KO_SMARCA2dTAG_dTAG47_24h_BRD4",
                    "SMARCA4KO_SMARCA2dTAG_DMSO_H3K27ac",
                    "SMARCA4KO_SMARCA2dTAG_dTAG47_24h_H3K27ac",
                    "SMARCC1KO_SMARCC2dTAG_DMSO_ARID1A",
                    "SMARCC1KO_SMARCC2dTAG_dTAG47_24h_ARID1A",
                    "SMARCC1KO_SMARCC2dTAG_DMSO_BRD4",
                    "SMARCC1KO_SMARCC2dTAG_dTAG47_24h_BRD4",
                    "SMARCC1KO_SMARCC2dTAG_DMSO_H3K27ac",
                    "SMARCC1KO_SMARCC2dTAG_dTAG47_24h_H3K27ac")
                    

data_macs2_plot <- data_macs2_qn_batch1[,columns_wanted]

rowsd_macs2_1<-apply(data_macs2_plot[,c(2:19)],1,sd)
rowmean_macs2_1<-rowMeans(data_macs2_plot[,c(2:19)])
Ztrans_macs2_1<-(data_macs2_plot[,c(2:19)]-rowmean_macs2_1)/rowsd_macs2_1

Ztrans_macs2_1$index <- data_macs2_plot$index

colnames(Ztrans_macs2_1)[1:18] <- paste(colnames(Ztrans_macs2_1)[1:18],"_batch1",sep = "")

### reading macs2 qn count data chip batch 2

data_macs2_qn_batch2 <- read.csv("Chip-Set2-quantile-norm-atac-dtag-consensus.csv")

colnames(data_macs2_qn_batch2)[1] <- "index"

columns_wanted <- c("index","WT_DMSO_H3K27ac","WT_SMARCA4dTAG_DMSO_H3K27ac","WT_SMARCA4dTAG_3h_dTAG_H3K27ac","WT_SMARCA4dTAG_6h_dTAG_H3K27ac","WT_SMARCA4dTAG_24h_dTAG_H3K27ac","WT_SMARCA4dTAG_72h_dTAG_H3K27ac","SMARCA4KO_DMSO_H3K27ac","SMARCA2KO_DMSO_H3K27ac","SMARCA4KO_SMARCA2dTAG_DMSO_H3K27ac","SMARCA4KO_SMARCA2dTAG_3h_dTAG_H3K27ac","SMARCA4KO_SMARCA2dTAG_6h_dTAG_H3K27ac","SMARCA4KO_SMARCA2dTAG_24h_dTAG_H3K27ac","SMARCA4KO_SMARCA2dTAG_72h_dTAG_H3K27ac")

data_macs2_plot_batch2 <- data_macs2_qn_batch2[,columns_wanted]

rowsd_macs2_2<-apply(data_macs2_plot_batch2[,c(2:14)],1,sd)
rowmean_macs2_2<-rowMeans(data_macs2_plot_batch2[,c(2:14)])
Ztrans_macs2_2<-(data_macs2_plot_batch2[,c(2:14)]-rowmean_macs2_2)/rowsd_macs2_2

Ztrans_macs2_2$index <- data_macs2_plot_batch2$index

colnames(Ztrans_macs2_2)[1:13] <- paste(colnames(Ztrans_macs2_2)[1:13],"_batch2",sep = "")


Ztrans_overall <- merge(Ztrans_macs2_1,Ztrans_macs2_2,"index")


### reading atac data

data <- read.table("BAF_Timecourse.matrix_cqn_norm_merged_replicates.tsv",header = T)

data_atac_want <- data[,c("WT_DMSO",                  
                          "WT_SMARCA4_F1_DMSO",
                          "WT_SMARCA4_F1_dTAG_47_2h",
                          "WT_SMARCA4_F1_dTAG_47_3h", 
                          "WT_SMARCA4_F1_dTAG_47_6h",
                          "WT_SMARCA4_F1_dTAG_47_24h",
                          "WT_SMARCA4_F1_dTAG_47_72h",
                          "SMARCA4_KO_DMSO",
                          "ARID2_KO_DMSO",
                          "A2_A4_F5_DMSO",
                          "A2_A4_F5_dTAG_47_2h",      
                          "A2_A4_F5_dTAG_47_3h",
                          "A2_A4_F5_dTAG_47_6h",
                          "A2_A4_F5_dTAG_47_24h",     
                          "A2_A4_F5_dTAG_47_72h",
                          "SMARCA2_KO_DMSO",          
                          "BRG_BRM_C8_DMSO",
                          "BRG_BRM_C8_dTAG_47_2h",
                          "BRG_BRM_C8_dTAG_47_3h",    
                          "BRG_BRM_C8_dTAG_47_6h",
                          "BRG_BRM_C8_dTAG_47_24h",
                          "BRG_BRM_C8_dTAG_47_72h",
                          "SMARCC2_KO_DMSO",
                          "SMARCC1_KO_DMSO",
                          "CC1_CC2_combined_control",
                          "CC1_CC2_F1_dTAG_47_3h",    
                          "CC1_CC2_F1_dTAG_47_6h","CC1_CC2_F1_dTAG_47_24h","CC1_CC2_F1_dTAG_47_72h")]



rowsd_atac<-apply(data_atac_want[,c(1:29)],1,sd)
rowmean_atac<-rowMeans(data_atac_want[,c(1:29)])
Ztrans_atac<-(data_atac_want[,c(1:29)]-rowmean_atac)/rowsd_atac

Ztrans_atac$index <- data$index


### reading in rna seq data

WT_SMARCA4dtag <- read.table("rna-seq-merged-vst-countsWT_SMARCA4dtag.txt",header = T)

WT_SMARCA4dtag_want <- WT_SMARCA4dtag[,c("gene_id","gene_name","HAP1_WT_DMSO_1_S67639_merged","HAP1_WT_SMARCA4_F1_DMSO_1_S67605_merged","HAP1_WT_SMARCA4_F1_dTAG_47_2h_1_S67612_merged","HAP1_WT_SMARCA4_F1_dTAG_47_3h_1_S67609_merged","HAP1_WT_SMARCA4_F1_dTAG_47_6h_1_S67610_merged","HAP1_WT_SMARCA4_F1_dTAG_47_24h_1_S67607_merged","HAP1_WT_SMARCA4_F1_dTAG_47_72h_1_S67608_merged","HAP1_SMARCA4_KO_DMSO_1_S67644_merged")]

ARID2KO_SMARCA4dtag <- read.table("rna-seq-merged-vst-countsARID2KO_SMARCA4dtag.txt",header = T)

ARID2KO_SMARCA4dtag_want <- ARID2KO_SMARCA4dtag[,c("gene_id","gene_name","HAP1_ARID2_KO_DMSO_1_S67642_merged","HAP1_A2_A4_F5_DMSO_1_S67600_merged","HAP1_A2_A4_F5_dTAG_47_2h_1_S67606_merged","HAP1_A2_A4_F5_dTAG_47_3h_1_S67603_merged","HAP1_A2_A4_F5_dTAG_47_6h_1_S67604_merged","HAP1_A2_A4_F5_dTAG_47_24h_1_S67601_merged","HAP1_A2_A4_F5_dTAG_47_72h_1_S67602_merged")]

BRG_BRMdtag <- read.table("rna-seq-merged-vst-countsBRG_BRMdtag.txt",header = T)

BRG_BRMdtag_want <- BRG_BRMdtag[,c("gene_id","gene_name","HAP1_SMARCA2_KO_DMSO_1_S67645_merged","HAP1_BRG_BRM_C8_DMSO_1_S67599_merged","HAP1_BRG_BRM_C8_dTAG_47_2h_1_S67638_merged","HAP1_BRG_BRM_C8_dTAG_47_3h_1_S67641_merged","HAP1_BRG_BRM_C8_dTAG_47_6h_1_S67640_merged","HAP1_BRG_BRM_C8_dTAG_47_24h_1_S67637_merged","HAP1_BRG_BRM_C8_dTAG_47_72h_1_S67636_merged")]

CC1_CC2dtag <- read.table("rna-seq-merged-vst-countsCC1_CC2dtag.txt",header = T)

CC1_CC2dtag_want <- CC1_CC2dtag[,c("gene_id","gene_name","HAP1_SMARCC2_KO_DMSO_1_S67648_merged","HAP1_SMARCC1_KO_DMSO_1_S67643_merged","HAP1_CC1_CC2_F1_DMSO_1_S67611_merged","HAP1_CC1_CC2_F1_dTAG_47_2h_1_S67598_merged","HAP1_CC1_CC2_F1_dTAG_47_3h_1_S67615_merged","HAP1_CC1_CC2_F1_dTAG_47_6h_1_S67616_merged","HAP1_CC1_CC2_F1_dTAG_47_24h_1_S67613_merged","HAP1_CC1_CC2_F1_dTAG_47_72h_1_S67614_merged")]

rna_list <- list(WT_SMARCA4dtag_want,ARID2KO_SMARCA4dtag_want,BRG_BRMdtag_want,CC1_CC2dtag_want)

rna_seq_counts <- Reduce(function(x,y) merge(x,y,"gene_id"),rna_list)

data_gene <- read.table("rna-atac-seq-merged-all-consensus-foldchanges-index.txt",header = T)

data_gene_want <- data_gene[,c("index","gene_name","gene_id")]

data_rna_counts <- merge(rna_seq_counts,data_gene_want,"gene_id")

rowsd_rna<-apply(data_rna_counts[,c(3:10,12:18,20:26,28:35)],1,sd)
rowmean_rna<-rowMeans(data_rna_counts[,c(3:10,12:18,20:26,28:35)])
Ztrans_rna<-(data_rna_counts[,c(3:10,12:18,20:26,28:35)]-rowmean_rna)/rowsd_rna

Ztrans_rna$index <- data_rna_counts$index

colnames(Ztrans_rna) <- gsub("HAP1_","",colnames(Ztrans_rna))

colnames(Ztrans_rna) <- gsub("_merged","",colnames(Ztrans_rna))

colnames(Ztrans_rna)[1:30] <- paste(colnames(Ztrans_rna)[1:30],"_rna",sep="")



### joining datasets

atac_chip <- merge(Ztrans_atac,Ztrans_overall,"index")

atac_chip_rna <- join(atac_chip,Ztrans_rna,"index")


### reading in cluster data FC based

cluster <- read.table("cluster_SL.annotated.bed")

colnames(cluster) <- c("index","cluster")

### getting data to plot dtag tc clusters

data_plot_cluster <- join(cluster,atac_chip_rna,"index")

### plotting heatmap

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

data_plot_heatmap <- data_plot_cluster[,-c(1,2)]

row.names(data_plot_heatmap) <- paste(data_plot_cluster$index,"_",row.names(data_plot_cluster),sep="")


ht2 <- Heatmap(data_plot_heatmap[,1:29], name = "Z-score counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 4),cluster_rows = F,cluster_columns = F,col=col_fun)

ht3 <- Heatmap(data_plot_heatmap[,30:47], name = "Z-score counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 4),cluster_rows = F,cluster_columns = F,col=col_fun)

ht4 <- Heatmap(data_plot_heatmap[,48:60], name = "Z-score counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 4),cluster_rows = F,cluster_columns = F,col=col_fun)

ht5 <- Heatmap(data_plot_heatmap[,61:90], name = "Z-score counts", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 4),cluster_rows = F,cluster_columns = F,col=col_fun)


pdf("SL-atac-dtag-clusters-w-counts-atac-chip-rna.pdf")
print(ht2 + ht3 + ht4 + ht5)
dev.off()
