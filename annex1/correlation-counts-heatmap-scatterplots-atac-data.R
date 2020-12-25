library(corrplot)
library(GGally)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(plyr)
library(RColorBrewer)
library(ggplot2)
library(hexbin)
library(scales)

### setting working directory

setwd("../data/")

count_data_merge <- read.table("BAF_Timecourse_Set3.matrix_cqn_norm_merged_replicates.tsv",header = T)

count_data_merge$test_dummy <- -(count_data_merge$A2_A4_F5_DMSO)

colnames(count_data_merge)[c(2:8,16,24:30)] <- c("HAP1_A2_A4_F5_DMSO_1_S67600_merged","HAP1_A2_A4_F5_dTAG_47_2h_1_S67606_merged",
                                "HAP1_A2_A4_F5_dTAG_47_3h_1_S67603_merged",
                                "HAP1_A2_A4_F5_dTAG_47_6h_1_S67604_merged","HAP1_A2_A4_F5_dTAG_47_24h_1_S67601_merged",
                                "HAP1_A2_A4_F5_dTAG_47_72h_1_S67602_merged","HAP1_ARID2_KO_DMSO_1_S67642_merged",
                                "HAP1_SMARCA4_KO_DMSO_1_S67644_merged","HAP1_WT_DMSO_1_S67639_merged.x",               
                                "HAP1_WT_SMARCA4_F1_DMSO_1_S67605_merged",
                                "HAP1_WT_SMARCA4_F1_dTAG_47_2h_1_S67612_merged","HAP1_WT_SMARCA4_F1_dTAG_47_3h_1_S67609_merged", 
                                "HAP1_WT_SMARCA4_F1_dTAG_47_6h_1_S67610_merged","HAP1_WT_SMARCA4_F1_dTAG_47_24h_1_S67607_merged","HAP1_WT_SMARCA4_F1_dTAG_47_72h_1_S67608_merged"
                                )

### getting order

correct_order <- c("HAP1_WT_DMSO_1_S67639_merged.x",               
                   "HAP1_WT_SMARCA4_F1_DMSO_1_S67605_merged",
                   "HAP1_WT_SMARCA4_F1_dTAG_47_2h_1_S67612_merged","HAP1_WT_SMARCA4_F1_dTAG_47_3h_1_S67609_merged", 
                  "HAP1_WT_SMARCA4_F1_dTAG_47_6h_1_S67610_merged","HAP1_WT_SMARCA4_F1_dTAG_47_24h_1_S67607_merged","HAP1_WT_SMARCA4_F1_dTAG_47_72h_1_S67608_merged",
                  "HAP1_SMARCA4_KO_DMSO_1_S67644_merged","HAP1_ARID2_KO_DMSO_1_S67642_merged","HAP1_A2_A4_F5_DMSO_1_S67600_merged",           
                  "HAP1_A2_A4_F5_dTAG_47_2h_1_S67606_merged",      
                    "HAP1_A2_A4_F5_dTAG_47_3h_1_S67603_merged","HAP1_A2_A4_F5_dTAG_47_6h_1_S67604_merged", "HAP1_A2_A4_F5_dTAG_47_24h_1_S67601_merged",     
                "HAP1_A2_A4_F5_dTAG_47_72h_1_S67602_merged", "test_dummy","index"  )


count_data_merge_interest <- count_data_merge[,correct_order]

colnames(count_data_merge_interest) <- gsub("_1_S67.*","",colnames(count_data_merge_interest))

colnames(count_data_merge_interest) <- gsub("HAP1_","",colnames(count_data_merge_interest))

colnames(count_data_merge_interest) <- gsub("_F1","",colnames(count_data_merge_interest))

colnames(count_data_merge_interest) <- gsub("A2_A4_F5","ARID2KO_SMARCA4",colnames(count_data_merge_interest))


### making correlation matrix

M_corr <- data.frame(cor(count_data_merge_interest[,-17],method = "pearson",use = "complete.obs"))

write.table(M_corr,"correlation-matrix-cqn-counts-atac-seq-all-genes-WT-SMARCA4dtag-and-ARID2KO-SMARCA4dtag-clones.txt",quote = F,sep = "\t")

matrix_corr <- as.matrix(M_corr)

col<- colorRampPalette(c("blue", "white", "red"))(10)

pdf("correlation-cqn-counts-atac-seq-all-genes-WT-SMARCA4dtag-and-ARID2KO-SMARCA4dtag-clones.pdf")
corrplot(matrix_corr,type = "lower",tl.cex = 0.5, col = col)
dev.off()

pdf("correlation-cqn-counts-atac-seq-all-genes-WT-SMARCA4dtag-and-ARID2KO-SMARCA4dtag-clones-R-value.pdf")
corrplot(matrix_corr,type = "lower",tl.cex = 0.5, method = "number",col = col,number.cex = 0.7)
dev.off()





### heatmap

count_data_merge_interest_wo_Inf_1 <- apply(count_data_merge_interest,2,function(x) gsub(Inf,NA,x))

count_data_merge_interest_wo_Inf_2 <- apply(count_data_merge_interest_wo_Inf_1,2,function(x) gsub(-Inf,NA,x))

count_data_merge_interest_numeric <- data.frame(apply(count_data_merge_interest_wo_Inf_2[,-17],2,function(x) as.numeric(x)))

count_data_merge_interest_numeric$index <- count_data_merge_interest_wo_Inf_2[,"index"]

data_cluster_regions <- read.table("cluster_all.annotated.bed",header = F)

colnames(data_cluster_regions) <- c("index","cluster")

count_data_merge_interest_numeric_all_cluster <- merge(data_cluster_regions,count_data_merge_interest_numeric,"index")

count_data_heatmap <- as.matrix(count_data_merge_interest_numeric_all_cluster[,-c(1,2,18)])

## plotting heatmap with Z-score

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

### calculating Zscore

rowsd<-apply(count_data_heatmap,1,sd)
rowmean<-rowMeans(count_data_heatmap)
Ztrans<-(count_data_heatmap-rowmean)/rowsd

## with own hclust function

d <- dist(Ztrans, method = "canberra") # distance matrix
fit <- hclust(d, method="ward.D")


ht1 <- Heatmap(Ztrans, name = "Z-score cqn counts", column_title = "Sample", row_title = "Gene", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = fit,cluster_columns = F,col=col_fun)

pdf("heatmap-atac-seq-Z-score-cqn-counts-clustered-canberra-ward-rows-all-11cluster-regions-ARID2KO-SMARCA4dtag-series-all-samples.pdf")
print(ht1)
dev.off()



ht2 <- Heatmap(Ztrans, name = "Z-score cqn counts", column_title = "Sample", row_title = "Gene", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = fit,col=col_fun)

pdf("heatmap-atac-seq-Z-score-cqn-counts-clustered-canberra-ward-rows-columns-all-11cluster-regions-ARID2KO-SMARCA4dtag-series-all-samples.pdf")
print(ht2)
dev.off()






## plotting heatmap without Z-score

col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))

## with own hclust function

d <- dist(count_data_heatmap, method = "canberra") # distance matrix
fit <- hclust(d, method="ward.D")


ht1 <- Heatmap(count_data_heatmap, name = "cqn counts", column_title = "Sample", row_title = "Gene", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = fit,cluster_columns = F,col=col_fun)

pdf("heatmap-atac-seq-cqn-counts-clustered-canberra-ward-rows-all-11cluster-regions-ARID2KO-SMARCA4dtag-series-all-samples.pdf")
print(ht1)
dev.off()



ht2 <- Heatmap(count_data_heatmap, name = "cqn counts", column_title = "Sample", row_title = "Gene", column_title_side = "bottom",row_names_gp = gpar(fontsize = 0.1), column_names_gp = gpar(fontsize = 6),cluster_rows = fit,col=col_fun)

pdf("heatmap-atac-seq-cqn-counts-clustered-canberra-ward-rows-columns-all-11cluster-regions-ARID2KO-SMARCA4dtag-series-all-samples.pdf")
print(ht2)
dev.off()


### scatter plots

cols = viridis(10)

scatter_data <- count_data_merge_interest[,-c(16,17)]

count1 = 1

count2 = 1

while (count1 < 16){

print(count1)  

while (count2 < 16){
  
print(count2)

correlation = cor(scatter_data[,count1],scatter_data[,count2],method = "pearson",use="na.or.complete")

scatter <- ggplot(scatter_data, aes(x=eval(parse(text = colnames(scatter_data[count1]))), y=eval(parse(text = colnames(scatter_data[count2]))))) + scale_fill_gradientn(colours = cols,limits=c(0,1000),oob=squish) + geom_hex(bins = 70) + theme(legend.position = "none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 


pdf(paste("scatterplot-atac-seq-all-11cluster-regions-",colnames(scatter_data[count1]),"-vs-",colnames(scatter_data[count2]),".pdf",sep=""))
print(scatter + annotate("text", x = 0, y = 10, label = c(paste("R = ",round(correlation,2),sep=""))) + xlab(colnames(scatter_data[count1])) + ylab(colnames(scatter_data[count2])) ) 
dev.off()


count2 <- count2 + 1

print(count2)
}
  
count1 <- count1 + 1

print(count1)

count2 <- count1 

print(count2)
   
}
