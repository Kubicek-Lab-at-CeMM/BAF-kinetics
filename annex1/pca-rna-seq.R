### getting packages

library(data.table)
library(factoextra)
library(NbClust)
library(stringr)
library(circlize)
library(colorspace)
library(GetoptLong)

### Setting working directory

setwd("../data/")

#### parental clones

data_vst <- read.table("rna-seq-merged-vst-countsWT_SMARCA4dtag.txt",header = T)

data_vst_ARID2KO <- read.table("rna-seq-merged-vst-countsARID2KO_SMARCA4dtag.txt",header = T)

data_vst_BRG_BRM <- read.table("rna-seq-merged-vst-countsBRG_BRMdtag.txt",header = T)

data_vst_CC1_CC2 <- read.table("rna-seq-merged-vst-countsCC1_CC2dtag.txt",header = T)

list_df <- list(data_vst,data_vst_ARID2KO,data_vst_BRG_BRM,data_vst_CC1_CC2)

data_merge <- Reduce(function(x,y) merge(x,y,"gene_id"),list_df)

data_merge_interest <- data_merge[,c("HAP1_WT_SMARCA4_F1_DMSO_1_S67605_merged","HAP1_SMARCA4_KO_DMSO_1_S67644_merged.x","HAP1_WT_DMSO_1_S67639_merged.x","HAP1_A2_A4_F5_DMSO_1_S67600_merged","HAP1_ARID2_KO_DMSO_1_S67642_merged","HAP1_WT_SMARCA4_F1_dTAG_47_24h_1_S67607_merged",
                                     "HAP1_WT_SMARCA4_F1_dTAG_47_2h_1_S67612_merged","HAP1_WT_SMARCA4_F1_dTAG_47_3h_1_S67609_merged",
                                     "HAP1_WT_SMARCA4_F1_dTAG_47_6h_1_S67610_merged","HAP1_WT_SMARCA4_F1_dTAG_47_72h_1_S67608_merged",
                                     "HAP1_A2_A4_F5_dTAG_47_24h_1_S67601_merged","HAP1_A2_A4_F5_dTAG_47_2h_1_S67606_merged",     
                                     "HAP1_A2_A4_F5_dTAG_47_3h_1_S67603_merged","HAP1_A2_A4_F5_dTAG_47_6h_1_S67604_merged",      
                                     "HAP1_A2_A4_F5_dTAG_47_72h_1_S67602_merged")]

data_merge_interest_clean <- data_merge_interest[! apply(data_merge_interest,1,function(x) sum(x)) == 0,]

colnames(data_merge_interest_clean) <- gsub("_1_S.*","",colnames(data_merge_interest_clean))



### PCA analysis

res.pca <- prcomp(t(data_merge_interest_clean), scale = TRUE)


pdf("pc1-2-rna-seq-WT-SMARCA4dtag-and-ARID2KO-SMARCA4dtag-series-clones-merged-rep.pdf",width = 7, height = 5)
fviz_pca_ind(res.pca,
             labelsize = 3,
             pointshape = 20,
             pointsize = 3,
             repel = T)
dev.off()



# separate series WT SMARCA4 dtag

data_merge_interest <- data_merge[,c("HAP1_WT_SMARCA4_F1_DMSO_1_S67605_merged","HAP1_SMARCA4_KO_DMSO_1_S67644_merged.x","HAP1_WT_DMSO_1_S67639_merged.x","HAP1_WT_SMARCA4_F1_dTAG_47_24h_1_S67607_merged",
                                     "HAP1_WT_SMARCA4_F1_dTAG_47_2h_1_S67612_merged","HAP1_WT_SMARCA4_F1_dTAG_47_3h_1_S67609_merged",
                                     "HAP1_WT_SMARCA4_F1_dTAG_47_6h_1_S67610_merged","HAP1_WT_SMARCA4_F1_dTAG_47_72h_1_S67608_merged"
                                    )]

data_merge_interest_clean <- data_merge_interest[! apply(data_merge_interest,1,function(x) sum(x)) == 0,]

colnames(data_merge_interest_clean) <- gsub("_S67.*","",colnames(data_merge_interest_clean))



### PCA analysis

res.pca <- prcomp(t(data_merge_interest_clean), scale = TRUE)


pdf("pc1-2-rna-seq-WT-SMARCA4dtag-series-clones-merged-rep.pdf",width = 7, height = 5)
fviz_pca_ind(res.pca,
             labelsize = 2,
             pointshape = 20,
             pointsize = 3,
             repel = T)
dev.off()


### merged replicates series ARID2KO SMARCA4dtag

data_merge_interest <- data_merge[,c("HAP1_SMARCA4_KO_DMSO_1_S67644_merged.x","HAP1_WT_DMSO_1_S67639_merged.x","HAP1_A2_A4_F5_DMSO_1_S67600_merged",
                                     "HAP1_ARID2_KO_DMSO_1_S67642_merged","HAP1_A2_A4_F5_dTAG_47_24h_1_S67601_merged","HAP1_A2_A4_F5_dTAG_47_2h_1_S67606_merged",     
                                     "HAP1_A2_A4_F5_dTAG_47_3h_1_S67603_merged","HAP1_A2_A4_F5_dTAG_47_6h_1_S67604_merged",      
                                     "HAP1_A2_A4_F5_dTAG_47_72h_1_S67602_merged")]

data_merge_interest_clean <- data_merge_interest[! apply(data_merge_interest,1,function(x) sum(x)) == 0,]

colnames(data_merge_interest_clean) <- gsub("_S67.*","",colnames(data_merge_interest_clean))



### PCA analysis

res.pca <- prcomp(t(data_merge_interest_clean), scale = TRUE)


pdf("pc1-2-rna-seq-ARID2KO-SMARCA4dtag-series-clones-merged-rep.pdf",width = 7, height = 5)
fviz_pca_ind(res.pca,
             labelsize = 2,
             pointshape = 20,
             pointsize = 3,
             repel = T)
dev.off()
