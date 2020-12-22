### getting packages

library(data.table)
library(factoextra)
library(NbClust)
library(stringr)
library(circlize)
library(colorspace)
library(GetoptLong)

### Setting working directory

setwd("../data/")

celllines <- c("ARID2KO_SMARCA4dtag","BRG_BRMdtag","CC1_CC2dtag","WT_SMARCA4dtag_C4","WT_SMARCA4dtag")

### reading in vst normalized count data

for (cellline in celllines){

print(cellline)
  
count = 7

data_vst <- read.table(paste("rnaseq_deseq_",cellline,"_counts_normalised.tsv",sep = ""),header = T) 

stopper = length(colnames(data_vst))

while (count < stopper){ 

print(count)
  
data_vst[,paste(colnames(data_vst)[count],"_merged",sep="")] <- (data_vst[,count] + data_vst[,count + 1]) / 2

count = count + 2

}

write.table(data_vst,paste("rna-seq-merged-vst-counts",cellline,".txt",sep = ""),quote = F,sep = "\t",row.names = F)

}


### reading differential genes

setwd("../data/")

genes_diff <- read.table("all-diff-genes.txt",header = T)

for (cellline in celllines){
  
  print(cellline)
 
  data <- read.table(paste("rna-seq-merged-vst-counts",cellline,".txt",sep = ""),header = T) 
  
  data_diff <- data[data$gene_id %in% genes_diff$gene_id,]
  
  data_diff_matrix <- as.matrix(data_diff[,str_detect(colnames(data_diff),"merged") == T])
  
  row.names(data_diff_matrix) <- data_diff$gene_id
  
  data_diff_matrix_clean <- data_diff_matrix[! apply(data_diff_matrix,1,function(x) sum(x)) == 0,]



#### PCA analysis

## WT SMARCA4dtag time-course

data_vst <- read.table("rna-seq-merged-vst-countsWT_SMARCA4dtag.txt",header = T)

data_merged <- data_vst[,str_detect(colnames(data_vst),"merged") == T]

data_clean <- data_merged[,-2]

data_diff_matrix_clean <- data_clean[! apply(data_clean,1,function(x) sum(x)) == 0,]

### PCA analysis

res.pca <- prcomp(t(data_diff_matrix_clean), scale = TRUE)

pdf(paste("pca-plot-PC1and2-WT_SMARCA4dtag.pdf",sep = ""),width = 7, height = 5,)
print(fviz_pca_ind(res.pca,
                   col.ind.sup = c("#5ac18e", "#0069b4", "#b51865","#6f1a16", "#00aee7", "#762a02", "#fcd116", "#518561"),
                   labelsize = 4,
                   pointshape = 20,
                   col.ind  = colnames(data_diff_matrix_clean),
                   repel = T,
                   label = F
))
dev.off()

