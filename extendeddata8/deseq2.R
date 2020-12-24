### Getting packages

library(data.table)
library(DESeq2)

setwd("../data/")

### ATP depletion time-course samples

### Reading raw count data

data <- read.csv("BAF-ATP-depletion.matrix_raw.csv")


### 6h BI sample

### making count table

countdata <- data[,c("R1_24h_DMSO", "R2_24h_DMSO","R1_5min_DMSO","R2_5min_DMSO","R1_6h_BI_protac","R2_6h_BI_protac")]                      
row.names(countdata) <- data$X

coldata <- data.frame(c("DMSO","DMSO","DMSO","DMSO","6h_BI","6h_BI"))
colnames(coldata) <- "condition"

row.names(coldata) <- c("R1_24h_DMSO", "R2_24h_DMSO","R1_5min_DMSO","R2_5min_DMSO","R1_6h_BI_protac","R2_6h_BI_protac")

### Generate DESeq object 

dds <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~ condition)


dds$condition <- relevel(dds$condition, ref = "DMSO")

analysis_res <- DESeq(dds)

results_res <- results(analysis_res)

results_res$peak <- data$X

results_res$comparison <- paste("6h_BI_vs_WT_DMSO",sep = "")

write.table(results_res, "Diff-Expression-6hBI-vs-WT-DMSO.txt",quote=F,sep="\t",row.names=F)

results_res_clean <- na.omit(results_res)

results_res_clean_sig <- results_res_clean[abs(results_res_clean$log2FoldChange) > 1 & results_res_clean$padj < 0.01,]

write.table(results_res_clean_sig,"Diff-Expression-6hBI-vs-WT-DMSO-sigregions.txt",quote=F,sep="\t",row.names=F)

### 24h BI sample

### making count table

countdata <- data[,c("R1_24h_DMSO", "R2_24h_DMSO","R1_5min_DMSO","R2_5min_DMSO","R1_24h_BI_protac","R2_24h_BI_protac")]                      
row.names(countdata) <- data$X

coldata <- data.frame(c("DMSO","DMSO","DMSO","DMSO","24h_BI","24h_BI"))
colnames(coldata) <- "condition"

row.names(coldata) <- c("R1_24h_DMSO", "R2_24h_DMSO","R1_5min_DMSO","R2_5min_DMSO","R1_24h_BI_protac","R2_24h_BI_protac")

### Generate DESeq object 

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ condition)


dds$condition <- relevel(dds$condition, ref = "DMSO")

analysis_res <- DESeq(dds)

results_res <- results(analysis_res)

results_res$peak <- data$X

results_res$comparison <- paste("24h_BI_vs_WT_DMSO",sep = "")

write.table(results_res, "Diff-Expression-24hBI-vs-WT-DMSO.txt",quote=F,sep="\t",row.names=F)


results_res_clean <- na.omit(results_res)

results_res_clean_sig <- results_res_clean[abs(results_res_clean$log2FoldChange) > 1 & results_res_clean$padj < 0.01,]

write.table(results_res_clean_sig,"Diff-Expression-24hBI-vs-WT-DMSO-sigregions.txt",quote=F,sep="\t",row.names=F)

