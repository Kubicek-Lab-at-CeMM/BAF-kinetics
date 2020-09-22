### getting packages

library(ggplot2)
library(hexbin)
library(viridis)
library(scales)

### setting working directory

setwd("../data/")

### reading in data

data_FC <- read.table("ATP-depletion-timecourse-allregions-new-deseq2.txt",header = T,stringsAsFactors = F)

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

list_samples <- list(Novartis_5min,Novartis_10min,Novartis_30min,Novartis_1h,Novartis_6h,Novartis_24h,BI_protac_6h,BI_protac_24h)

data_ATP <- Reduce(function(x,y) merge(x,y,"index"),list_samples)

colnames(data_ATP)[1] <- "index_compound"

### reading in SMARCA4 data

SMARCA4_data <- read.csv("differential_analysis.deseq_result.all_comparisons.csv")

WT_SMARCA4dtag_2h <- SMARCA4_data[SMARCA4_data$comparison_name == "ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_2h_vs_DMSO",]

colnames(WT_SMARCA4dtag_2h)[2:8] <- paste(colnames(WT_SMARCA4dtag_2h)[2:8],"_SMARCA4dtag_2h",sep="")

WT_SMARCA4dtag_3h <- SMARCA4_data[SMARCA4_data$comparison_name == "ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_3h_vs_DMSO",]

colnames(WT_SMARCA4dtag_3h)[2:8] <- paste(colnames(WT_SMARCA4dtag_3h)[2:8],"_SMARCA4dtag_3h",sep="")

WT_SMARCA4dtag_6h <- SMARCA4_data[SMARCA4_data$comparison_name == "ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_6h_vs_DMSO",]

colnames(WT_SMARCA4dtag_6h)[2:8] <- paste(colnames(WT_SMARCA4dtag_6h)[2:8],"_SMARCA4dtag_6h",sep="")

WT_SMARCA4dtag_24h <- SMARCA4_data[SMARCA4_data$comparison_name == "ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_24h_vs_DMSO",]

colnames(WT_SMARCA4dtag_24h)[2:8] <- paste(colnames(WT_SMARCA4dtag_24h)[2:8],"_SMARCA4dtag_24h",sep="")

WT_SMARCA4dtag_72h <- SMARCA4_data[SMARCA4_data$comparison_name == "ATAC-seq_HAP1_WT_SMARCA4_F1_dTAG_47_72h_vs_DMSO",]

colnames(WT_SMARCA4dtag_72h)[2:8] <- paste(colnames(WT_SMARCA4dtag_72h)[2:8],"_SMARCA4dtag_72h",sep="")

WT_SMARCA4dtag_KO <- SMARCA4_data[SMARCA4_data$comparison_name == "ATAC-seq_HAP1_BRGKO_vs_DMSO",]

colnames(WT_SMARCA4dtag_KO)[2:8] <- paste(colnames(WT_SMARCA4dtag_KO)[2:8],"_SMARCA4dtag_KO",sep="")

### merging with ATP depletion data

list_samples <- list(WT_SMARCA4dtag_2h,WT_SMARCA4dtag_3h,WT_SMARCA4dtag_6h,WT_SMARCA4dtag_24h,WT_SMARCA4dtag_72h,WT_SMARCA4dtag_KO)

data_SMARCA4 <- Reduce(function(x,y) merge(x,y,"index"),list_samples)
  
colnames(data_SMARCA4)[1] <- "index_dtag"

## getting merge file

merge_file <- read.table("dtag_ATP_overlap_merge.bed",header=F,stringsAsFactors = F)

colnames(merge_file) <- c("chr_dtag","pos1_dtag","pos2_dtag","chr_compound","pos1_compound","pos2_compound")

merge_file$index_dtag <- paste(merge_file$chr_dtag,":",merge_file$pos1_dtag,"-",merge_file$pos2_dtag,sep="")

merge_file$index_compound <- paste(merge_file$chr_compound,":",merge_file$pos1_compound,"-",merge_file$pos2_compound,sep="")

data_merge <- merge(data_ATP,merge_file,"index_compound")

data_merge2 <- merge(data_merge,data_SMARCA4,"index_dtag")



### plotting correlation for all regions 

cols = viridis(10)

### Novartis compounds

correlation = cor(data_merge2$log2FoldChange_N_24h,data_merge2$log2FoldChange_SMARCA4dtag_24h,method = "pearson",use="na.or.complete")

scatter <- ggplot(data_merge2, aes(x=log2FoldChange_N_24h, y=log2FoldChange_SMARCA4dtag_24h)) + scale_fill_gradientn(colours = cols,limits=c(0,1000),oob=squish) + geom_hex(bins = 70) + ylim(-9,9) + xlim(-9,9) + theme(legend.position = "none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 

pdf(paste("scatterplot-all-regions-SMARCA4dtag24h-Novartis24h.pdf",sep=""))
print(scatter + annotate("text", x = -6, y = 6, label = c(paste("R = ",round(correlation,2),sep="")))) 
dev.off()
