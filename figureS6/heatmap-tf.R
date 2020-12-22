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

cluster_vector <- c("cluster1","cluster2","cluster3","cluster4","cluster5")

for (element in cluster_vector){

### reading in data

homer_known_list[[element]] <- fread(paste("",element,"/knownResults.txt",sep = ""),data.table = F)

homer_known_list[[element]]$ID <- paste(homer_known_list[[element]]$`Motif Name`,"_",homer_known_list[[element]]$Consensus,sep = "")

colnames(homer_known_list[[element]])[c(1:9)] <- paste(colnames(homer_known_list[[element]])[c(1:9)],"_",element,sep = "")

}

### merging data frames

homer_known_merge <- Reduce(function(x, y) merge(x, y, by = "ID", all=TRUE), homer_known_list)

#homer_known_filtered <- homer_known_merge[homer_known_merge$`P-value_cluster1` < 1e-10 | homer_known_merge$`P-value_cluster2` < 1e-10 | homer_known_merge$`P-value_cluster3` < 1e-10 | homer_known_merge$`P-value_cluster4` < 1e-10 | homer_known_merge$`P-value_cluster5` < 1e-10,]

#write.table(homer_known_filtered,"200714-homer-results-known-ATP-depletion-regions-filteredp10.txt",quote = F, sep = "\t",row.names = F)

homer_known_filtered <- homer_known_merge[homer_known_merge$`P-value_cluster1` < 1e-20 | homer_known_merge$`P-value_cluster2` < 1e-20 | homer_known_merge$`P-value_cluster3` < 1e-20 | homer_known_merge$`P-value_cluster4` < 1e-20 | homer_known_merge$`P-value_cluster5` < 1e-20,]

#write.table(homer_known_filtered,"200714-homer-results-known-ATP-depletion-regions-filteredp20.txt",quote = F, sep = "\t",row.names = F)

#homer_known_filtered <- homer_known_merge[homer_known_merge$`P-value_cluster1` < 1e-25 | homer_known_merge$`P-value_cluster2` < 1e-25 | homer_known_merge$`P-value_cluster3` < 1e-25 | homer_known_merge$`P-value_cluster4` < 1e-25 | homer_known_merge$`P-value_cluster5` < 1e-25,]

#write.table(homer_known_filtered,"200714-homer-results-known-ATP-depletion-regions-filteredp25.txt",quote = F, sep = "\t",row.names = F)


#homer_known_filtered <- homer_known_merge[homer_known_merge$`P-value_cluster1` < 1e-40 | homer_known_merge$`P-value_cluster2` < 1e-40 | homer_known_merge$`P-value_cluster3` < 1e-40 | homer_known_merge$`P-value_cluster4` < 1e-40 | homer_known_merge$`P-value_cluster5` < 1e-40,]

#write.table(homer_known_filtered,"200714-homer-results-known-known-ATP-depletion-regions-filteredp40.txt",quote = F, sep = "\t",row.names = F)


### making plot table from p values

plot_table <- homer_known_filtered[,c("P-value_cluster1","P-value_cluster2","P-value_cluster3","P-value_cluster4","P-value_cluster5")]

row.names(plot_table) <- homer_known_filtered$`Motif Name_cluster1`

plot_table_log <- apply(plot_table,2,function(x) -log10(x))

### plotting heatmap

cols <- brewer.pal(9, "Oranges")

pal <- colorRampPalette(cols)

col_fun = colorRamp2(c(0,15,30,45,60,75),pal(6))

ht_just_p <- Heatmap(plot_table_log, name = "-log10 p", column_title = "Cluster", row_title = "Feature", column_title_side = "bottom",row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 14),cluster_columns = F,clustering_distance_rows = "euclidean",col=col_fun)

pdf("Heatmap-homer-ATP-depl-clusters-filtered-p20.pdf")
print(ht_just_p)
dev.off()

