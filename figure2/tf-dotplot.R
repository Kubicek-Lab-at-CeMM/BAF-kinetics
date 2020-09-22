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

homer_known_filtered <- homer_known_merge[homer_known_merge$`P-value_cluster1_tc` < 1e-10 | homer_known_merge$`P-value_cluster2_tc` < 1e-10 | homer_known_merge$`P-value_cluster3_tc` < 1e-10 | homer_known_merge$`P-value_cluster1_cluster2_KO` < 1e-10 | homer_known_merge$`P-value_cluster3_KO` < 1e-10,]

write.table(homer_known_filtered,"homer-results-known-5-cluster-filteredp10.txt",quote = F, sep = "\t",row.names = F)

homer_known_filtered <- homer_known_merge[homer_known_merge$`P-value_cluster1_tc` < 1e-20 | homer_known_merge$`P-value_cluster2_tc` < 1e-20 | homer_known_merge$`P-value_cluster3_tc` < 1e-20 | homer_known_merge$`P-value_cluster1_cluster2_KO` < 1e-20 | homer_known_merge$`P-value_cluster3_KO` < 1e-20,]

write.table(homer_known_filtered,"homer-results-known-5-cluster-filteredp20.txt",quote = F, sep = "\t",row.names = F)

homer_known_filtered <- homer_known_merge[homer_known_merge$`P-value_cluster1_tc` < 1e-40 | homer_known_merge$`P-value_cluster2_tc` < 1e-40 | homer_known_merge$`P-value_cluster3_tc` < 1e-40 | homer_known_merge$`P-value_cluster1_cluster2_KO` < 1e-40 | homer_known_merge$`P-value_cluster3_KO` < 1e-40,]

write.table(homer_known_filtered,"homer-results-known-5-cluster-filteredp40.txt",quote = F, sep = "\t",row.names = F)


### making plot table from p values

plot_table <- homer_known_filtered[,c("P-value_cluster1_tc","P-value_cluster2_tc","P-value_cluster3_tc","P-value_cluster1_cluster2_KO","P-value_cluster3_KO","% of Target Sequences with Motif_cluster1_tc","% of Target Sequences with Motif_cluster2_tc","% of Target Sequences with Motif_cluster3_tc","% of Target Sequences with Motif_cluster1_cluster2_KO","% of Target Sequences with Motif_cluster3_KO")]

row.names(plot_table) <- homer_known_filtered$`Motif Name_cluster1_tc`

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

plot_table_stacked$cluster[1:18] <- "cluster1"

plot_table_stacked$cluster[19:36] <- "cluster2"

plot_table_stacked$cluster[37:54] <- "cluster3"

plot_table_stacked$cluster[55:72] <- "cluster4"

plot_table_stacked$cluster[73:90] <- "cluster5"

plot_table_stacked$feature <- rep(row.names(plot_table_log),5)


cols <- c("deepskyblue",
          "red")


dotplot = ggplot(plot_table_stacked) +
  geom_point(mapping = aes(x= cluster, y=feature, size=Motif_Ratio, colour = `-log10(p_value)`)) +
  scale_colour_gradientn(colours = cols, breaks=c(0,40,80), limits = c(0,80),oob=squish,na.value = 'black') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

pdf("dotplot-homer-known-SMARCA4-timecourse-KO-clusters-p20.pdf",useDingbats = F,width = 10, height = 5) 
print(dotplot)
dev.off()
