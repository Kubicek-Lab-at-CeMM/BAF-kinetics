library(data.table)
library(factoextra)
library(NbClust)
library(circlize)
library(colorspace)
library(GetoptLong)

### setting working directory

setwd("../data/")

### reading data

data_rpm_batch3 <- read.csv("BAF_Timecourse_Set3.matrix_norm_rpm.csv")

data_rpm_want <- data_rpm_batch3[,c("ATAC.seq_HAP1_A2_A4_F5_DMSO_1",            
                   "ATAC.seq_HAP1_A2_A4_F5_dTAG_47_24h_1","ATAC.seq_HAP1_A2_A4_F5_dTAG_47_2h_1",      
                   "ATAC.seq_HAP1_A2_A4_F5_dTAG_47_3h_1","ATAC.seq_HAP1_A2_A4_F5_dTAG_47_6h_1",      
                   "ATAC.seq_HAP1_A2_A4_F5_dTAG_47_72h_1","ATAC.seq_HAP1_A2_A4_F5_DMSO_2",            
                   "ATAC.seq_HAP1_A2_A4_F5_dTAG_47_24h_2",      "ATAC.seq_HAP1_A2_A4_F5_dTAG_47_2h_2",      
                   "ATAC.seq_HAP1_A2_A4_F5_dTAG_47_3h_2",       "ATAC.seq_HAP1_A2_A4_F5_dTAG_47_6h_2",      
                   "ATAC.seq_HAP1_A2_A4_F5_dTAG_47_72h_2",      "ATAC.seq_HAP1_ARID2_KO_DMSO_1","ATAC.seq_HAP1_ARID2_KO_DMSO_2",
                   "ATAC.seq_HAP1_SMARCA4_KO_DMSO_1",          
                   "ATAC.seq_HAP1_SMARCA4_KO_DMSO_2",
                   "ATAC.seq_HAP1_WT_SMARCA4_F1_DMSO_1",        "ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_1",
                   "ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_1",  "ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_1", 
                   "ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_1",  "ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_1",
                   "ATAC.seq_HAP1_WT_SMARCA4_F1_DMSO_2",       "ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_2",
                   "ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_2",  "ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_2", 
                   "ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_2",  "ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_2",
                   "ATAC.seq_HAP1_WT_DMSO_1",                   "ATAC.seq_HAP1_WT_DMSO_2")]

### merging counts

data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_DMSO <- (data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_DMSO_1 + data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_DMSO_2) / 2

data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_dTAG_47_24h <- (data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_dTAG_47_24h_1 + data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_dTAG_47_24h_2) / 2

data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_dTAG_47_2h <- (data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_dTAG_47_2h_1 + data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_dTAG_47_2h_2) / 2

data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_dTAG_47_72h <- (data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_dTAG_47_72h_1 + data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_dTAG_47_72h_2) / 2

data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_dTAG_47_6h <- (data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_dTAG_47_6h_1 + data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_dTAG_47_6h_2) / 2

data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_dTAG_47_3h <- (data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_dTAG_47_3h_1 + data_rpm_want$ATAC.seq_HAP1_A2_A4_F5_dTAG_47_3h_2) / 2

data_rpm_want$ATAC.seq_HAP1_ARID2_KO_DMSO <- (data_rpm_want$ATAC.seq_HAP1_ARID2_KO_DMSO_1 + data_rpm_want$ATAC.seq_HAP1_ARID2_KO_DMSO_2) / 2

data_rpm_want$ATAC.seq_HAP1_SMARCA4_KO_DMSO <- (data_rpm_want$ATAC.seq_HAP1_SMARCA4_KO_DMSO_1 + data_rpm_want$ATAC.seq_HAP1_SMARCA4_KO_DMSO_2) / 2

data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_DMSO <- (data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_DMSO_1 + data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_DMSO_2) / 2

data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h <- (data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_1 + data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_2) / 2

data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h <- (data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_1 + data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_2) / 2

data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h <- (data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_1 + data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_2) / 2

data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h <- (data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_1 + data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_2) / 2

data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h <- (data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_1 + data_rpm_want$ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_2) / 2

data_rpm_want$ATAC.seq_HAP1_WT_DMSO <- (data_rpm_want$ATAC.seq_HAP1_WT_DMSO_1 + data_rpm_want$ATAC.seq_HAP1_WT_DMSO_2) / 2


### PCA analysis

res.pca <- prcomp(t(data_rpm_want[,31:45]), scale = TRUE)


pdf("pc1-2-atac-seq-WT-SMARCA4dtag-and-ARID2KO-SMARCA4dtag-series-clones-merged-rep.pdf",width = 7, height = 5)
fviz_pca_ind(res.pca,
             labelsize = 2,
             pointshape = 20,
             pointsize = 3,
             repel = T)
dev.off()



### WT SMARCA4 dtag series

data_rpm_want_1 <- data_rpm_want[,c(
                                    "ATAC.seq_HAP1_SMARCA4_KO_DMSO",
                                    "ATAC.seq_HAP1_WT_SMARCA4_F1_DMSO",        "ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h",
                                    "ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h",  "ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h", 
                                    "ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h",  "ATAC.seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h",
                                    "ATAC.seq_HAP1_WT_DMSO")]

### PCA analysis

res.pca <- prcomp(t(data_rpm_want_1), scale = TRUE)


pdf("pc1-2-atac-seq-WT-SMARCA4dtag-series-clones-merged-rep.pdf",width = 7, height = 5)
fviz_pca_ind(res.pca,
             labelsize = 2,
             pointshape = 20,
             pointsize = 3,
             repel = T)
dev.off()



# ARID2KO SMARCA4dtag

data_rpm_want_2 <- data_rpm_want[,c("ATAC.seq_HAP1_A2_A4_F5_DMSO",            
                                    "ATAC.seq_HAP1_A2_A4_F5_dTAG_47_24h","ATAC.seq_HAP1_A2_A4_F5_dTAG_47_2h",      
                                    "ATAC.seq_HAP1_A2_A4_F5_dTAG_47_3h","ATAC.seq_HAP1_A2_A4_F5_dTAG_47_6h",      
                                    "ATAC.seq_HAP1_A2_A4_F5_dTAG_47_72h","ATAC.seq_HAP1_ARID2_KO_DMSO",
                                    "ATAC.seq_HAP1_SMARCA4_KO_DMSO",          
                          "ATAC.seq_HAP1_WT_DMSO")]

### PCA analysis

res.pca <- prcomp(t(data_rpm_want_2), scale = TRUE)


pdf("pc1-2-atac-seq-ARID2KO-SMARCA4dtag-series-clones-merged-rep.pdf",width = 7, height = 5)
fviz_pca_ind(res.pca,
             labelsize = 2,
             pointshape = 20,
             pointsize = 3,
             repel = T)
dev.off()


