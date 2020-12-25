## getting packages

library(data.table)
library(DESeq2)

## setting working directory

setwd("/Users/sgrosche/Development/201031-atac-seq-analyses/DESeq2/")

## getting file

data <- read.csv("../data/BAF_Timecourse_Set3.matrix_raw.txt")

### WT DMSO vs ARID2KO SMARCA4dtag DMSO

### making count table

countdata <- data[,c("ATAC.seq_HAP1_A2_A4_F5_DMSO_1", "ATAC.seq_HAP1_A2_A4_F5_DMSO_2","ATAC.seq_HAP1_WT_DMSO_1","ATAC.seq_HAP1_WT_DMSO_2")]       

row.names(countdata) <- paste(data$X,sep="")

coldata <- data.frame(c("A2_A4_F5_DMSO","A2_A4_F5_DMSO","WT_DMSO","WT_DMSO"))
colnames(coldata) <- "condition"

row.names(coldata) <- c("ATAC.seq_HAP1_A2_A4_F5_DMSO_1", "ATAC.seq_HAP1_A2_A4_F5_DMSO_2","ATAC.seq_HAP1_WT_DMSO_1","ATAC.seq_HAP1_WT_DMSO_2")

### Generate DESeq object 

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ condition)


dds$condition <- relevel(dds$condition, ref = "A2_A4_F5_DMSO")

analysis_res <- DESeq(dds)

results_res <- results(analysis_res)

results_res$peak <- paste(data$X,sep="")

results_res$comparison <- "WT_DMSO_vs_A2_A4_F5_DMSO"

write.table(results_res, paste("201101-diff-expression-WT_DMSO_vs_A2_A4_F5_DMSO.txt",sep=""),quote=F,sep="\t",row.names=F)



results_res_clean <- na.omit(results_res)

results_res_clean_sig <- results_res_clean[abs(results_res_clean$log2FoldChange) > 1 & results_res_clean$padj < 0.01,]

write.table(results_res_clean_sig,paste("diff-expression-WT_DMSO_vs_A2_A4_F5_DMSO-sig-regions.txt",sep = ""),quote=F,sep="\t",row.names=F)


# volcano

### Reading in data

data_res <- data.frame(results_res)

data_res$neg_log10_p <- -log10(data_res$padj)

data_res$coloring <- ifelse(data_res$padj < 0.01 & abs(data_res$log2FoldChange) > 1,'sig',"not sig")


volcano_plot <- ggplot(data_res) + geom_point(aes(x= log2FoldChange, y=neg_log10_p, color = coloring)) + 
  scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() + ylim(0,40)

pdf(paste("volcano-plot-ATAC-WT-DMSO-ARID2KO-SMARCA4dtag-DMSO-all-consensus-scaled.pdf"))
print(volcano_plot)
dev.off()

### getting up and down significant regions

data_diff_sig <- na.omit(data_res[data_res$coloring == "sig", ])

write.table(data_diff_sig,"sig-regions-WT-DMSO-ARID2KO-SMARCA4dtag-DMSO.txt",quote = F,sep = "\t",row.names = F)

data_diff_up <- na.omit(data_res[data_res$coloring == "sig" & data_res$log2FoldChange > 1, ])

write.table(data_diff_up,"sig-upregulated-regions-WT-DMSO-ARID2KO-SMARCA4dtag-DMSO.txt",quote = F,sep = "\t",row.names = F)

data_diff_down <- na.omit(data_res[data_res$coloring == "sig" & data_res$log2FoldChange < -1, ])

write.table(data_diff_down,"sig-downregulated-regions-WT-DMSO-ARID2KO-SMARCA4dtag-DMSO.txt",quote = F,sep = "\t",row.names = F)


### WT SMARCA4dtag vs SMARCA4KO

### making count table

countdata <- data[,c("ATAC.seq_HAP1_SMARCA4_KO_DMSO_1", "ATAC.seq_HAP1_SMARCA4_KO_DMSO_2","ATAC.seq_HAP1_WT_SMARCA4_F1_DMSO_1","ATAC.seq_HAP1_WT_SMARCA4_F1_DMSO_2")]       

row.names(countdata) <- paste(data$X,sep="")

coldata <- data.frame(c("SMARCA4_KO_DMSO","SMARCA4_KO_DMSO","WT_SMARCA4_F1_DMSO","WT_SMARCA4_F1_DMSO"))
colnames(coldata) <- "condition"

row.names(coldata) <- c("ATAC.seq_HAP1_SMARCA4_KO_DMSO_1", "ATAC.seq_HAP1_SMARCA4_KO_DMSO_2","ATAC.seq_HAP1_WT_SMARCA4_F1_DMSO_1","ATAC.seq_HAP1_WT_SMARCA4_F1_DMSO_2")

### Generate DESeq object 

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ condition)


dds$condition <- relevel(dds$condition, ref = "SMARCA4_KO_DMSO")

analysis_res <- DESeq(dds)

results_res <- results(analysis_res)

results_res$peak <- paste(data$X,sep="")

results_res$comparison <- "WT_SMARCA4dtag_DMSO_vs_SMARCA4_KO_DMSO"

write.table(results_res, paste("diff-expression-WT_SMARCA4dtag_DMSO_vs_SMARCA4_KO_DMSO.txt",sep=""),quote=F,sep="\t",row.names=F)



results_res_clean <- na.omit(results_res)

results_res_clean_sig <- results_res_clean[abs(results_res_clean$log2FoldChange) > 1 & results_res_clean$padj < 0.01,]

write.table(results_res_clean_sig,paste("diff-expression-WT_SMARCA4dtag_DMSO_vs_SMARCA4_KO_DMSO-sig-regions.txt",sep = ""),quote=F,sep="\t",row.names=F)


# volcano

### Reading in data

data_res <- data.frame(results_res)

data_res$neg_log10_p <- -log10(data_res$padj)

data_res$coloring <- ifelse(data_res$padj < 0.01 & abs(data_res$log2FoldChange) > 1,'sig',"not sig")


volcano_plot <- ggplot(data_res) + geom_point(aes(x= log2FoldChange, y=neg_log10_p, color = coloring)) + 
  scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() + ylim(0,40)

pdf(paste("volcano-plot-ATAC-WT_SMARCA4dtag_DMSO_vs_SMARCA4_KO_DMSO-all-consensus-scaled.pdf"))
print(volcano_plot)
dev.off()

### getting up and down significant regions

data_diff_sig <- na.omit(data_res[data_res$coloring == "sig", ])

write.table(data_diff_sig,"sig-regions-WT_SMARCA4dtag_DMSO_vs_SMARCA4_KO_DMSO.txt",quote = F,sep = "\t",row.names = F)

data_diff_up <- na.omit(data_res[data_res$coloring == "sig" & data_res$log2FoldChange > 1, ])

write.table(data_diff_up,"sig-upregulated-regions-WT_SMARCA4dtag_DMSO_vs_SMARCA4_KO_DMSO.txt",quote = F,sep = "\t",row.names = F)

data_diff_down <- na.omit(data_res[data_res$coloring == "sig" & data_res$log2FoldChange < -1, ])

write.table(data_diff_down,"sig-downregulated-regions-WT_SMARCA4dtag_DMSO_vs_SMARCA4_KO_DMSO.txt",quote = F,sep = "\t",row.names = F)



### WT vs SMARCA4KO

### making count table

countdata <- data[,c("ATAC.seq_HAP1_SMARCA4_KO_DMSO_1", "ATAC.seq_HAP1_SMARCA4_KO_DMSO_2","ATAC.seq_HAP1_WT_DMSO_1","ATAC.seq_HAP1_WT_DMSO_2")]       

row.names(countdata) <- paste(data$X,sep="")

coldata <- data.frame(c("SMARCA4_KO_DMSO","SMARCA4_KO_DMSO","WT_DMSO","WT_DMSO"))
colnames(coldata) <- "condition"

row.names(coldata) <- c("ATAC.seq_HAP1_SMARCA4_KO_DMSO_1", "ATAC.seq_HAP1_SMARCA4_KO_DMSO_2","ATAC.seq_HAP1_WT_DMSO_1","ATAC.seq_HAP1_WT_DMSO_2")

### Generate DESeq object 

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ condition)


dds$condition <- relevel(dds$condition, ref = "SMARCA4_KO_DMSO")

analysis_res <- DESeq(dds)

results_res <- results(analysis_res)

results_res$peak <- paste(data$X,sep="")

results_res$comparison <- "WT_DMSO_vs_SMARCA4_KO_DMSO"

write.table(results_res, paste("diff-expression-WT_DMSO_vs_SMARCA4_KO_DMSO.txt",sep=""),quote=F,sep="\t",row.names=F)



results_res_clean <- na.omit(results_res)

results_res_clean_sig <- results_res_clean[abs(results_res_clean$log2FoldChange) > 1 & results_res_clean$padj < 0.01,]

write.table(results_res_clean_sig,paste("diff-expression-WT_DMSO_vs_SMARCA4_KO_DMSO-sig-regions.txt",sep = ""),quote=F,sep="\t",row.names=F)


# volcano

### Reading in data

data_res <- data.frame(results_res)

data_res$neg_log10_p <- -log10(data_res$padj)

data_res$coloring <- ifelse(data_res$padj < 0.01 & abs(data_res$log2FoldChange) > 1,'sig',"not sig")


volcano_plot <- ggplot(data_res) + geom_point(aes(x= log2FoldChange, y=neg_log10_p, color = coloring)) + 
  scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() + ylim(0,40)

pdf(paste("volcano-plot-ATAC-WT_DMSO_vs_SMARCA4_KO_DMSO-all-consensus-scaled.pdf"))
print(volcano_plot)
dev.off()

### getting up and down significant regions

data_diff_sig <- na.omit(data_res[data_res$coloring == "sig", ])

write.table(data_diff_sig,"sig-regions-WT_DMSO_vs_SMARCA4_KO_DMSO.txt",quote = F,sep = "\t",row.names = F)

data_diff_up <- na.omit(data_res[data_res$coloring == "sig" & data_res$log2FoldChange > 1, ])

write.table(data_diff_up,"sig-upregulated-regions-WT_DMSO_vs_SMARCA4_KO_DMSO.txt",quote = F,sep = "\t",row.names = F)

data_diff_down <- na.omit(data_res[data_res$coloring == "sig" & data_res$log2FoldChange < -1, ])

write.table(data_diff_down,"sig-downregulated-regions-WT_DMSO_vs_SMARCA4_KO_DMSO.txt",quote = F,sep = "\t",row.names = F)














### SMARCA4KO vs WT

### making count table

countdata <- data[,c("ATAC.seq_HAP1_SMARCA4_KO_DMSO_1", "ATAC.seq_HAP1_SMARCA4_KO_DMSO_2","ATAC.seq_HAP1_WT_DMSO_1","ATAC.seq_HAP1_WT_DMSO_2")]       

row.names(countdata) <- paste(data$X,sep="")

coldata <- data.frame(c("SMARCA4_KO_DMSO","SMARCA4_KO_DMSO","WT_DMSO","WT_DMSO"))
colnames(coldata) <- "condition"

row.names(coldata) <- c("ATAC.seq_HAP1_SMARCA4_KO_DMSO_1", "ATAC.seq_HAP1_SMARCA4_KO_DMSO_2","ATAC.seq_HAP1_WT_DMSO_1","ATAC.seq_HAP1_WT_DMSO_2")

### Generate DESeq object 

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ condition)


dds$condition <- relevel(dds$condition, ref = "WT_DMSO")

analysis_res <- DESeq(dds)

results_res <- results(analysis_res)

results_res$peak <- paste(data$X,sep="")

results_res$comparison <- "SMARCA4_KO_DMSO_vs_WT_DMSO"

write.table(results_res, paste("201101-diff-expression-SMARCA4_KO_DMSO_vs_WT_DMSO.txt",sep=""),quote=F,sep="\t",row.names=F)



results_res_clean <- na.omit(results_res)

results_res_clean_sig <- results_res_clean[abs(results_res_clean$log2FoldChange) > 1 & results_res_clean$padj < 0.01,]

write.table(results_res_clean_sig,paste("201101-diff-expression-SMARCA4_KO_DMSO_vs_WT_DMSO-sig-regions.txt",sep = ""),quote=F,sep="\t",row.names=F)


# volcano

### Reading in data

data_res <- data.frame(results_res)

data_res$neg_log10_p <- -log10(data_res$padj)

data_res$coloring <- ifelse(data_res$padj < 0.01 & abs(data_res$log2FoldChange) > 1,'sig',"not sig")


volcano_plot <- ggplot(data_res) + geom_point(aes(x= log2FoldChange, y=neg_log10_p, color = coloring)) + 
  scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() + ylim(0,40) + xlim(-10,5)

pdf(paste("201101-volcano-plot-ATAC-SMARCA4_KO_DMSO_vs_WT_DMSO-all-consensus-scaled.pdf"))
print(volcano_plot)
dev.off()

### getting up and down significant regions

data_diff_sig <- na.omit(data_res[data_res$coloring == "sig", ])

write.table(data_diff_sig,"201101-sig-regions-SMARCA4_KO_DMSO_vs_WT_DMSO.txt",quote = F,sep = "\t",row.names = F)

data_diff_up <- na.omit(data_res[data_res$coloring == "sig" & data_res$log2FoldChange > 1, ])

write.table(data_diff_up,"201101-sig-upregulated-regions-SMARCA4_KO_DMSO_vs_WT_DMSO.txt",quote = F,sep = "\t",row.names = F)

data_diff_down <- na.omit(data_res[data_res$coloring == "sig" & data_res$log2FoldChange < -1, ])

write.table(data_diff_down,"201101-sig-downregulated-regions-SMARCA4_KO_DMSO_vs_WT_DMSO.txt",quote = F,sep = "\t",row.names = F)









### ARID2KO SMARCA4dtag vs SMARCA4KO

### making count table

countdata <- data[,c("ATAC.seq_HAP1_SMARCA4_KO_DMSO_1", "ATAC.seq_HAP1_SMARCA4_KO_DMSO_2","ATAC.seq_HAP1_A2_A4_F5_DMSO_1","ATAC.seq_HAP1_A2_A4_F5_DMSO_2")]       

row.names(countdata) <- paste(data$X,sep="")

coldata <- data.frame(c("SMARCA4_KO_DMSO","SMARCA4_KO_DMSO","ARID2KO_SMARCA4dtag","ARID2KO_SMARCA4dtag"))
colnames(coldata) <- "condition"

row.names(coldata) <- c("ATAC.seq_HAP1_SMARCA4_KO_DMSO_1", "ATAC.seq_HAP1_SMARCA4_KO_DMSO_2","ATAC.seq_HAP1_A2_A4_F5_DMSO_1","ATAC.seq_HAP1_A2_A4_F5_DMSO_2")

### Generate DESeq object 

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ condition)


dds$condition <- relevel(dds$condition, ref = "SMARCA4_KO_DMSO")

analysis_res <- DESeq(dds)

results_res <- results(analysis_res)

results_res$peak <- paste(data$X,sep="")

results_res$comparison <- "ARID2KO_SMARCA4dtag_vs_SMARCA4_KO_DMSO"

write.table(results_res, paste("201101-diff-expression-ARID2KO_SMARCA4dtag_vs_SMARCA4_KO_DMSO.txt",sep=""),quote=F,sep="\t",row.names=F)



results_res_clean <- na.omit(results_res)

results_res_clean_sig <- results_res_clean[abs(results_res_clean$log2FoldChange) > 1 & results_res_clean$padj < 0.01,]

write.table(results_res_clean_sig,paste("201101-diff-expression-ARID2KO_SMARCA4dtag_vs_SMARCA4_KO_DMSO-sig-regions.txt",sep = ""),quote=F,sep="\t",row.names=F)


# volcano

### Reading in data

data_res <- data.frame(results_res)

data_res$neg_log10_p <- -log10(data_res$padj)

data_res$coloring <- ifelse(data_res$padj < 0.01 & abs(data_res$log2FoldChange) > 1,'sig',"not sig")


volcano_plot <- ggplot(data_res) + geom_point(aes(x= log2FoldChange, y=neg_log10_p, color = coloring)) + 
  scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() + ylim(0,40)

pdf(paste("201101-volcano-plot-ATAC-ARID2KO_SMARCA4dtag_vs_SMARCA4_KO_DMSO-all-consensus-scaled.pdf"))
print(volcano_plot)
dev.off()

### getting up and down significant regions

data_diff_sig <- na.omit(data_res[data_res$coloring == "sig", ])

write.table(data_diff_sig,"201101-sig-regions-ARID2KO_SMARCA4dtag_vs_SMARCA4_KO_DMSO.txt",quote = F,sep = "\t",row.names = F)

data_diff_up <- na.omit(data_res[data_res$coloring == "sig" & data_res$log2FoldChange > 1, ])

write.table(data_diff_up,"201101-sig-upregulated-regions-ARID2KO_SMARCA4dtag_vs_SMARCA4_KO_DMSO.txt",quote = F,sep = "\t",row.names = F)

data_diff_down <- na.omit(data_res[data_res$coloring == "sig" & data_res$log2FoldChange < -1, ])

write.table(data_diff_down,"201101-sig-downregulated-regions-ARID2KO_SMARCA4dtag_vs_SMARCA4_KO_DMSO.txt",quote = F,sep = "\t",row.names = F)







### ARID2KO vs SMARCA4KO

### making count table

countdata <- data[,c("ATAC.seq_HAP1_SMARCA4_KO_DMSO_1", "ATAC.seq_HAP1_SMARCA4_KO_DMSO_2","ATAC.seq_HAP1_ARID2_KO_DMSO_1","ATAC.seq_HAP1_ARID2_KO_DMSO_1")]       

row.names(countdata) <- paste(data$X,sep="")

coldata <- data.frame(c("SMARCA4_KO_DMSO","SMARCA4_KO_DMSO","ARID2KO","ARID2KO"))
colnames(coldata) <- "condition"

row.names(coldata) <- c("ATAC.seq_HAP1_SMARCA4_KO_DMSO_1", "ATAC.seq_HAP1_SMARCA4_KO_DMSO_2","ATAC.seq_HAP1_ARID2_KO_DMSO_1","ATAC.seq_HAP1_ARID2_KO_DMSO_2")

### Generate DESeq object 

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ condition)


dds$condition <- relevel(dds$condition, ref = "SMARCA4_KO_DMSO")

analysis_res <- DESeq(dds)

results_res <- results(analysis_res)

results_res$peak <- paste(data$X,sep="")

results_res$comparison <- "ARID2KO_vs_SMARCA4_KO_DMSO"

write.table(results_res, paste("201101-diff-expression-ARID2KO_vs_SMARCA4_KO_DMSO.txt",sep=""),quote=F,sep="\t",row.names=F)



results_res_clean <- na.omit(results_res)

results_res_clean_sig <- results_res_clean[abs(results_res_clean$log2FoldChange) > 1 & results_res_clean$padj < 0.01,]

write.table(results_res_clean_sig,paste("201101-diff-expression-ARID2KO_vs_SMARCA4_KO_DMSO-sig-regions.txt",sep = ""),quote=F,sep="\t",row.names=F)


# volcano

### Reading in data

data_res <- data.frame(results_res)

data_res$neg_log10_p <- -log10(data_res$padj)

data_res$coloring <- ifelse(data_res$padj < 0.01 & abs(data_res$log2FoldChange) > 1,'sig',"not sig")


volcano_plot <- ggplot(data_res) + geom_point(aes(x= log2FoldChange, y=neg_log10_p, color = coloring)) + 
  scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() + ylim(0,40)

pdf(paste("201101-volcano-plot-ATAC-ARID2KO_vs_SMARCA4_KO_DMSO-all-consensus-scaled.pdf"))
print(volcano_plot)
dev.off()

### getting up and down significant regions

data_diff_sig <- na.omit(data_res[data_res$coloring == "sig", ])

write.table(data_diff_sig,"201101-sig-regions-ARID2KO_vs_SMARCA4_KO_DMSO.txt",quote = F,sep = "\t",row.names = F)

data_diff_up <- na.omit(data_res[data_res$coloring == "sig" & data_res$log2FoldChange > 1, ])

write.table(data_diff_up,"201101-sig-upregulated-regions-ARID2KO_vs_SMARCA4_KO_DMSO.txt",quote = F,sep = "\t",row.names = F)

data_diff_down <- na.omit(data_res[data_res$coloring == "sig" & data_res$log2FoldChange < -1, ])

write.table(data_diff_down,"201101-sig-downregulated-regions-ARID2KO_vs_SMARCA4_KO_DMSO.txt",quote = F,sep = "\t",row.names = F)



### SMARCA4dtag DMSO vs WT DMSO

### making count table

countdata <- data[,c("ATAC.seq_HAP1_WT_SMARCA4_F1_DMSO_1", "ATAC.seq_HAP1_WT_SMARCA4_F1_DMSO_2","ATAC.seq_HAP1_WT_DMSO_1","ATAC.seq_HAP1_WT_DMSO_2")]       

row.names(countdata) <- paste(data$X,sep="")

coldata <- data.frame(c("WT_SMARCA4_F1_DMSO","WT_SMARCA4_F1_DMSO","WT_DMSO","WT_DMSO"))
colnames(coldata) <- "condition"

row.names(coldata) <- c("ATAC.seq_HAP1_WT_SMARCA4_F1_DMSO_1", "ATAC.seq_HAP1_WT_SMARCA4_F1_DMSO_2","ATAC.seq_HAP1_WT_DMSO_1","ATAC.seq_HAP1_WT_DMSO_2")

### Generate DESeq object 

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ condition)


dds$condition <- relevel(dds$condition, ref = "WT_DMSO")

analysis_res <- DESeq(dds)

results_res <- results(analysis_res)

results_res$peak <- paste(data$X,sep="")

results_res$comparison <- "WT_SMARCA4_dtag_vs_WT_DMSO"

write.table(results_res, paste("201101-diff-expression-WT_SMARCA4_dtag_vs_WT_DMSO.txt",sep=""),quote=F,sep="\t",row.names=F)



results_res_clean <- na.omit(results_res)

results_res_clean_sig <- results_res_clean[abs(results_res_clean$log2FoldChange) > 1 & results_res_clean$padj < 0.01,]

write.table(results_res_clean_sig,paste("201101-diff-expression-WT_SMARCA4_dtag_vs_WT_DMSO-sig-regions.txt",sep = ""),quote=F,sep="\t",row.names=F)


# volcano

### Reading in data

data_res <- data.frame(results_res)

data_res$neg_log10_p <- -log10(data_res$padj)

data_res$coloring <- ifelse(data_res$padj < 0.01 & abs(data_res$log2FoldChange) > 1,'sig',"not sig")


volcano_plot <- ggplot(data_res) + geom_point(aes(x= log2FoldChange, y=neg_log10_p, color = coloring)) + 
  scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() + ylim(0,40)

pdf(paste("201101-volcano-plot-ATAC-WT_SMARCA4_dtag_vs_WT_DMSO-all-consensus-scaled.pdf"))
print(volcano_plot)
dev.off()

### getting up and down significant regions

data_diff_sig <- na.omit(data_res[data_res$coloring == "sig", ])

write.table(data_diff_sig,"201101-sig-regions-WT_SMARCA4_dtag_vs_WT_DMSO.txt",quote = F,sep = "\t",row.names = F)

data_diff_up <- na.omit(data_res[data_res$coloring == "sig" & data_res$log2FoldChange > 1, ])

write.table(data_diff_up,"201101-sig-upregulated-regions-WT_SMARCA4_dtag_vs_WT_DMSO.txt",quote = F,sep = "\t",row.names = F)

data_diff_down <- na.omit(data_res[data_res$coloring == "sig" & data_res$log2FoldChange < -1, ])

write.table(data_diff_down,"201101-sig-downregulated-regions-WT_SMARCA4_dtag_vs_WT_DMSO.txt",quote = F,sep = "\t",row.names = F)




### ARID2KO SMARCA4dtag vs ARID2KO

### making count table

countdata <- data[,c("ATAC.seq_HAP1_A2_A4_F5_DMSO_1", "ATAC.seq_HAP1_A2_A4_F5_DMSO_2","ATAC.seq_HAP1_ARID2_KO_DMSO_1","ATAC.seq_HAP1_ARID2_KO_DMSO_1")]       

row.names(countdata) <- paste(data$X,sep="")

coldata <- data.frame(c("ARID2KO_SMARCA4dtag","ARID2KO_SMARCA4dtag","ARID2KO","ARID2KO"))
colnames(coldata) <- "condition"

row.names(coldata) <- c("ATAC.seq_HAP1_A2_A4_F5_DMSO_1", "ATAC.seq_HAP1_A2_A4_F5_DMSO_2","ATAC.seq_HAP1_ARID2_KO_DMSO_1","ATAC.seq_HAP1_ARID2_KO_DMSO_1")

### Generate DESeq object 

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ condition)


dds$condition <- relevel(dds$condition, ref = "ARID2KO")

analysis_res <- DESeq(dds)

results_res <- results(analysis_res)

results_res$peak <- paste(data$X,sep="")

results_res$comparison <- "ARID2KO_SMARCA4dtag_vs_ARID2KO"

write.table(results_res, paste("201101-diff-expression-ARID2KO_SMARCA4dtag_vs_ARID2KO.txt",sep=""),quote=F,sep="\t",row.names=F)



results_res_clean <- na.omit(results_res)

results_res_clean_sig <- results_res_clean[abs(results_res_clean$log2FoldChange) > 1 & results_res_clean$padj < 0.01,]

write.table(results_res_clean_sig,paste("201101-diff-expression-ARID2KO_SMARCA4dtag_vs_ARID2KO-sig-regions.txt",sep = ""),quote=F,sep="\t",row.names=F)


# volcano

### Reading in data

data_res <- data.frame(results_res)

data_res$neg_log10_p <- -log10(data_res$padj)

data_res$coloring <- ifelse(data_res$padj < 0.01 & abs(data_res$log2FoldChange) > 1,'sig',"not sig")


volcano_plot <- ggplot(data_res) + geom_point(aes(x= log2FoldChange, y=neg_log10_p, color = coloring)) + 
  scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() + ylim(0,40)

pdf(paste("201101-volcano-plot-ATAC-ARID2KO_SMARCA4dtag_vs_ARID2KO-all-consensus-scaled.pdf"))
print(volcano_plot)
dev.off()

### getting up and down significant regions

data_diff_sig <- na.omit(data_res[data_res$coloring == "sig", ])

write.table(data_diff_sig,"201101-sig-regions-ARID2KO_SMARCA4dtag_vs_ARID2KO.txt",quote = F,sep = "\t",row.names = F)

data_diff_up <- na.omit(data_res[data_res$coloring == "sig" & data_res$log2FoldChange > 1, ])

write.table(data_diff_up,"201101-sig-upregulated-regions-ARID2KO_SMARCA4dtag_vs_ARID2KO.txt",quote = F,sep = "\t",row.names = F)

data_diff_down <- na.omit(data_res[data_res$coloring == "sig" & data_res$log2FoldChange < -1, ])

write.table(data_diff_down,"201101-sig-downregulated-regions-ARID2KO_SMARCA4dtag_vs_ARID2KO.txt",quote = F,sep = "\t",row.names = F)


