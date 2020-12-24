### installing packages

library(data.table)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(RColorBrewer)

### setting working directory

setwd("../data/")

### creating list 

homer_known_list <- list()

### looping over clusters

cluster_vector <- c("cluster1_tc","cluster2_tc","cluster3_tc","cluster1_cluster2_KO","cluster3_KO")

for (element in cluster_vector){

### reading in data

homer_known_list[[element]] <- fread(paste("",element,"/knownResults.txt",sep = ""),data.table = F)

homer_known_list[[element]]$ID <- paste(homer_known_list[[element]]$`Motif Name`,"_",homer_known_list[[element]]$Consensus,sep = "")

colnames(homer_known_list[[element]])[c(1:9)] <- paste(colnames(homer_known_list[[element]])[c(1:9)],"_",element,sep = "")

}

### merging data frames

homer_known_merge <- Reduce(function(x, y) merge(x, y, by = "ID", all=TRUE), homer_known_list)

motifs_to_plot <- homer_known_merge[unique(c(31,252,250,255,361,   364,    162,362,    358,  396)),]

### making plot table from p values

plot_table <- motifs_to_plot[,c("P-value_cluster1_tc","P-value_cluster2_tc","P-value_cluster3_tc","P-value_cluster1_cluster2_KO","P-value_cluster3_KO","% of Target Sequences with Motif_cluster1_tc","% of Target Sequences with Motif_cluster2_tc","% of Target Sequences with Motif_cluster3_tc","% of Target Sequences with Motif_cluster1_cluster2_KO","% of Target Sequences with Motif_cluster3_KO")]

row.names(plot_table) <- motifs_to_plot$`Motif Name_cluster1_tc`

plot_table_log <- plot_table

plot_table_log[1:5] <- apply(plot_table[1:5],2,function(x) -log10(x))

plot_table_log[6:10] <- apply(plot_table[6:10],2,function(x) gsub("%","",x))

plot_table_log[6:10] <- apply(plot_table_log[6:10],2,function(x) as.numeric(x))

plot_table_log[6:10] <- apply(plot_table_log[6:10],2,function(x) (x)/100)

cluster1 <- plot_table_log[,c("P-value_cluster1_tc","% of Target Sequences with Motif_cluster1_tc")]

cluster2 <- plot_table_log[,c("P-value_cluster2_tc","% of Target Sequences with Motif_cluster2_tc")]

cluster3 <- plot_table_log[,c("P-value_cluster3_tc","% of Target Sequences with Motif_cluster3_tc")]

cluster4 <- plot_table_log[,c("P-value_cluster1_cluster2_KO","% of Target Sequences with Motif_cluster1_cluster2_KO")]

cluster5 <- plot_table_log[,c("P-value_cluster3_KO","% of Target Sequences with Motif_cluster3_KO")]


df_list <- list(cluster1,cluster2,cluster3,cluster4,cluster5)

plot_table_stacked <- rbindlist(df_list,use.names = F)

colnames(plot_table_stacked) <- c("-log10(p_value)","Motif_Ratio")

plot_table_stacked$cluster <- NA

plot_table_stacked$cluster[1:10] <- "cluster1"

plot_table_stacked$cluster[11:20] <- "cluster2"

plot_table_stacked$cluster[21:30] <- "cluster3"

plot_table_stacked$cluster[31:40] <- "cluster4"

plot_table_stacked$cluster[41:50] <- "cluster5"

plot_table_stacked$feature <- rep(row.names(plot_table_log),5)

order_vector <- plot_table_stacked[order(plot_table_stacked$`-log10(p_value)`,decreasing = T),'feature']

plot_table_stacked$feature <- factor(plot_table_stacked$feature, levels = rev(unique(order_vector$feature)))


cols <- c("deepskyblue",
          "red")


dotplot = ggplot(plot_table_stacked) +
  geom_point(mapping = aes(x= cluster, y=feature, size=Motif_Ratio, colour = `-log10(p_value)`)) +
  scale_colour_gradientn(colours = cols, breaks=c(0,40,80), limits = c(0,80),oob=squish,na.value = 'black') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

pdf("dotplot-homer-candidates.pdf",useDingbats = F,width = 10, height = 5) 
print(dotplot)
dev.off()



