### installing packages

library(ggplot2)

### Setting working directory

setwd("../data/")

### 6h BI protac

### Reading in data

data_diff <- read.table("Diff-Expression-6hBI-vs-WT-DMSO.txt",header = T)

data_diff$neg_log10_p <- -log10(data_diff$padj)

data_diff$coloring <- ifelse(data_diff$padj < 0.01 & abs(data_diff$log2FoldChange) > 1,'sig',"not sig")


volcano_plot <- ggplot(data_diff) + geom_point(aes(x= log2FoldChange, y=neg_log10_p, color = coloring)) + 
   scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() 

pdf("volcano-plot-6h-BI-protac-vs-WT.pdf")
print(volcano_plot)
dev.off()


### 24h BI protac

### Reading in data

data_diff <- read.table("Diff-Expression-24hBI-vs-WT-DMSO.txt",header = T)

data_diff$neg_log10_p <- -log10(data_diff$padj)

data_diff$coloring <- ifelse(data_diff$padj < 0.01 & abs(data_diff$log2FoldChange) > 1,'sig',"not sig")


volcano_plot <- ggplot(data_diff) + geom_point(aes(x= log2FoldChange, y=neg_log10_p, color = coloring)) + 
  scale_color_manual(values=c("#000000", "#CC0000")) + theme(legend.title = element_blank()) + theme_classic() 

pdf("volcano-plot-24h-BI-protac-vs-WT.pdf")
print(volcano_plot)
dev.off()

