### installing packages

library(reshape2)
library(ggplot2)
library(scales)

### setting working directory

setwd("../data/")

### plots of timecourse data against all consensus

clustervector <- c("cluster1","cluster2","cluster3","cluster4","cluster5")

### making lists

list_df_data <- list()

### getting significant GO terms

sig_GO <- c("GO:0072148","GO:0003279","GO:0001708","GO:0007389","GO:0045165","GO:0000978","GO:0000987","GO:0001227","GO:0005096","GO:0072001","GO:0001822","GO:0001655")


### starting loop

for (element in clustervector){

### reading in data

cluster = element

print(cluster)

data = read.table(paste("",cluster,"/go-enrichment-results-",cluster,"-against-consensus.txt",sep = ""),header = T)


list_df_data[[cluster]] = data[data$ID %in% sig_GO,]

}

### making plot table

data_merge <- rbind(list_df_data[['cluster1']],list_df_data[['cluster4']])

data_merge$cluster <- NA

data_merge$cluster[1:9] <- "cluster1"

data_merge$cluster[10:12] <- "cluster4"

data_merge$gene_ratio_calc <- as.numeric(unlist(strsplit(as.character(data_merge$GeneRatio),"/"))[c(T,F)]) / as.numeric(unlist(strsplit(as.character(data_merge$GeneRatio),"/"))[c(F,T)])

data_merge$neg_log10_padj <- -log10(data_merge$p.adjust)

### making plot

cols <- c("red",
          "deepskyblue")


dotplot = ggplot(data_merge) +
  geom_point(mapping = aes(x= cluster, y=Description, size=gene_ratio_calc, colour = p.adjust)) +
  scale_colour_gradientn(colours = cols, breaks= c(0.01,0.02,0.03,0.04,0.05), limits = c(0,0.05),oob=squish) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) 

pdf("dotplot-GO-enrichment-timecourse-clusters-against-all-consensus-deepskyblue.pdf",useDingbats = F) 
print(dotplot)
dev.off() 
