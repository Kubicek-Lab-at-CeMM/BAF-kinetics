### installing packages

library(corrplot)
library(data.table)
library(factoextra)
library(NbClust)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)

### setting working directory

setwd("../data/")

### reading in data

data_rpm <- read.csv2("BAF-ATP-depletion.matrix_norm.rpm.csv",sep=",",stringsAsFactors = F)

colnames(data_rpm)[1] = "index"

data_rpm[2:length(colnames(data_rpm))] <- apply(data_rpm[2:length(colnames(data_rpm))],2,function(x) as.numeric(x))

### merging replicates rpm data

data_rpm$merge_24h_DMSO <- (data_rpm$R1_24h_DMSO + data_rpm$R2_24h_DMSO) /2

data_rpm$merge_24h_N <- (data_rpm$R1_24h_N + data_rpm$R2_24h_N) /2

data_rpm$merge_6h_N <- (data_rpm$R1_6h_N + data_rpm$R2_6h_N) /2

data_rpm$merge_1h_N <- (data_rpm$R1_1h_N + data_rpm$R2_1h_N) /2

data_rpm$merge_30min_N <- (data_rpm$R1_30min_N + data_rpm$R2_30min_N) /2

data_rpm$merge_10min_N <- (data_rpm$R1_10min_N + data_rpm$R2_10min_N) /2

data_rpm$merge_5min_N <- (data_rpm$R1_5min_N + data_rpm$R2_5min_N) /2

data_rpm$merge_5min_DMSO <- (data_rpm$R1_5min_DMSO + data_rpm$R2_5min_DMSO) /2

data_rpm$merge_24h_control <- (data_rpm$R1_24h_control + data_rpm$R2_24h_control) /2

data_rpm$merge_24h_BI_protac <- (data_rpm$R1_24h_BI_protac + data_rpm$R2_24h_BI_protac) /2

data_rpm$merge_6h_BI_protac <- (data_rpm$R1_6h_BI_protac + data_rpm$R2_6h_BI_protac) /2

data_rpm$merge_6h_control <- (data_rpm$R1_6h_control + data_rpm$R2_6h_control) /2

data_rpm$merge_combined_N_control <- (data_rpm$R1_24h_DMSO + data_rpm$R2_24h_DMSO + data_rpm$R1_5min_DMSO + data_rpm$R2_5min_DMSO) /4

write.table(data_rpm[,c(1,26:38)],"BAF-ATP-depletion.matrix_rpm_merged.txt",quote = F,sep = "\t",row.names = F)

### principle component analyses

res.pca <- prcomp(t(data_rpm[,c(26:33,38)]), scale = TRUE)


pdf("pca-plot-PC1and2-merged-replicates-rpm-counts-Novartis-inhibitor.pdf")
fviz_pca_ind(res.pca,
             labelsize = 4,
             pointshape = 20,
             col.ind = rownames(res.pca$x),
             repel = T)
dev.off()
