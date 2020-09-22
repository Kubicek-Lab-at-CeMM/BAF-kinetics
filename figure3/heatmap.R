### Installing packages

library(data.table)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(plyr)

### setting working directory

setwd("../data/")

### reading in data

data_FC <- read.table("ATP-depletion-timecourse-all-diff-regions-new-deseq2.txt",header = T,stringsAsFactors = F)

### data wrangling

Novartis_10min <- data_FC[data_FC$comparison_name == "10min_Novartis_inhibitor",]

colnames(Novartis_10min)[2:8] <- paste(colnames(Novartis_10min)[2:8],"_N_10min",sep="")

Novartis_5min <- data_FC[data_FC$comparison_name == "5min_Novartis_inhibitor",]

colnames(Novartis_5min)[2:8] <- paste(colnames(Novartis_5min)[2:8],"_N_5min",sep="")

Novartis_30min <- data_FC[data_FC$comparison_name == "30min_Novartis_inhibitor",]

colnames(Novartis_30min)[2:8] <- paste(colnames(Novartis_30min)[2:8],"_N_30min",sep="")

Novartis_1h <- data_FC[data_FC$comparison_name == "1h_Novartis_inhibitor",]

colnames(Novartis_1h)[2:8] <- paste(colnames(Novartis_1h)[2:8],"_N_1h",sep="")

Novartis_6h <- data_FC[data_FC$comparison_name == "6h_Novartis_inhibitor",]

colnames(Novartis_6h)[2:8] <- paste(colnames(Novartis_6h)[2:8],"_N_6h",sep="")

Novartis_24h <- data_FC[data_FC$comparison_name == "24h_Novartis_inhibitor",]

colnames(Novartis_24h)[2:8] <- paste(colnames(Novartis_24h)[2:8],"_N_24h",sep="")

BI_protac_6h <- data_FC[data_FC$comparison_name == "6h_BI_vs_WT_DMSO",]

colnames(BI_protac_6h)[2:8] <- paste(colnames(BI_protac_6h)[2:8],"_BI_protac_6h",sep="")

BI_protac_24h <- data_FC[data_FC$comparison_name == "24h_BI_vs_WT_DMSO",]

colnames(BI_protac_24h)[2:8] <- paste(colnames(BI_protac_24h)[2:8],"_BI_protac_24h",sep="")



### merging data

data_ATP <- cbind(Novartis_5min,Novartis_10min,Novartis_30min,Novartis_1h,Novartis_6h,Novartis_24h,BI_protac_6h,BI_protac_24h)

colnames(data_ATP)[1] <- "index_compound"

write.table(data_ATP,"ATP-depletion-data-new-deseq2-treatment-per-column-diff-regions.txt",quote = F,sep = "\t",row.names=F)

### chosing columns

data_plot <- as.matrix(data_ATP[,c(3,11,19,27,35,43,51,59)])

row.names(data_plot) <- data_ATP$index_compound

## plotting heatmap

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))


## with own hclust function
d <- dist(data_plot[,c(1:6)], method = "canberra") # distance matrix
fit <- hclust(d, method="ward.D")


ht2 <- Heatmap(data_plot[,c(1:6)], name = "log2 FC", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = fit,split = 5,cluster_columns = F,col=col_fun)

row_ordering <- unlist(row_order(ht2))

ht3 <- Heatmap(data_plot[,c(7:8)], name = "log2 FC", column_title = "Sample", row_title = "Genomic region", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = F,row_order = row_ordering,cluster_columns = F,col=col_fun)

pdf("Heatmap-FC-clustered-Novartis-inhibitor-tc.pdf")
print(ht2 + ht3)
dev.off()




### getting regions per cluster

cluster1 <- data.frame(row.names(data_plot[unlist(row_order(ht2)[1]),]))

colnames(cluster1) <- "index"

cluster1$chr <- unlist(strsplit(as.character(cluster1$index),":"))[c(T,F)]

cluster1$position <- unlist(strsplit(as.character(cluster1$index),":"))[c(F,T)]

cluster1$pos1 <- unlist(strsplit(as.character(cluster1$position),"-"))[c(T,F)]

cluster1$pos2 <- unlist(strsplit(as.character(cluster1$position),"-"))[c(F,T)]

cluster1_bed <- cluster1[,c("chr","pos1","pos2")]

write.table(cluster1_bed,"cluster1-FC-clustered.bed",quote=F,sep = "\t",row.names = F, col.names = F)

cluster2 <- data.frame(row.names(data_plot[unlist(row_order(ht2)[2]),]))

colnames(cluster2) <- "index"

cluster2$chr <- unlist(strsplit(as.character(cluster2$index),":"))[c(T,F)]

cluster2$position <- unlist(strsplit(as.character(cluster2$index),":"))[c(F,T)]

cluster2$pos1 <- unlist(strsplit(as.character(cluster2$position),"-"))[c(T,F)]

cluster2$pos2 <- unlist(strsplit(as.character(cluster2$position),"-"))[c(F,T)]

cluster2_bed <- cluster2[,c("chr","pos1","pos2")]

write.table(cluster2_bed,"cluster2-FC-clustered.bed",quote=F,sep = "\t",row.names = F, col.names = F)


cluster3 <- data.frame(row.names(data_plot[unlist(row_order(ht2)[3]),]))

colnames(cluster3) <- "index"

cluster3$chr <- unlist(strsplit(as.character(cluster3$index),":"))[c(T,F)]

cluster3$position <- unlist(strsplit(as.character(cluster3$index),":"))[c(F,T)]

cluster3$pos1 <- unlist(strsplit(as.character(cluster3$position),"-"))[c(T,F)]

cluster3$pos2 <- unlist(strsplit(as.character(cluster3$position),"-"))[c(F,T)]

cluster3_bed <- cluster3[,c("chr","pos1","pos2")]

write.table(cluster3_bed,"cluster3-FC-clustered.bed",quote=F,sep = "\t",row.names = F, col.names = F)


cluster4 <- data.frame(row.names(data_plot[unlist(row_order(ht2)[4]),]))

colnames(cluster4) <- "index"

cluster4$chr <- unlist(strsplit(as.character(cluster4$index),":"))[c(T,F)]

cluster4$position <- unlist(strsplit(as.character(cluster4$index),":"))[c(F,T)]

cluster4$pos1 <- unlist(strsplit(as.character(cluster4$position),"-"))[c(T,F)]

cluster4$pos2 <- unlist(strsplit(as.character(cluster4$position),"-"))[c(F,T)]

cluster4_bed <- cluster4[,c("chr","pos1","pos2")]

write.table(cluster4_bed,"cluster4-FC-clustered.bed",quote=F,sep = "\t",row.names = F, col.names = F)



cluster5 <- data.frame(row.names(data_plot[unlist(row_order(ht2)[5]),]))

colnames(cluster5) <- "index"

cluster5$chr <- unlist(strsplit(as.character(cluster5$index),":"))[c(T,F)]

cluster5$position <- unlist(strsplit(as.character(cluster5$index),":"))[c(F,T)]

cluster5$pos1 <- unlist(strsplit(as.character(cluster5$position),"-"))[c(T,F)]

cluster5$pos2 <- unlist(strsplit(as.character(cluster5$position),"-"))[c(F,T)]

cluster5_bed <- cluster5[,c("chr","pos1","pos2")]

write.table(cluster5_bed,"cluster5-FC-clustered.bed",quote=F,sep = "\t",row.names = F, col.names = F)



