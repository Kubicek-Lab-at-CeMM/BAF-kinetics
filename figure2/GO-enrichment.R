### Installing packages

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

 BiocManager::install("clusterProfiler")

 if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 
 BiocManager::install("DOSE")
 
 library(clusterProfiler)

 library(enrichplot)

 library(org.Hs.eg.db)

 library(DOSE)


### setting working directory
 
setwd("../data/")
 
### Go enrichment

for (element in (1:5)){

genes <- read.table(paste("genes-cluster",element,".txt",sep=""),header = T)

genes_clean <- genes[!genes$gene == "nan" & !genes$gene == "" & !genes$gene == ".","gene"]

universe <- read.table("all-genes-consensus.txt",header = T)

universe_clean <- universe[!universe$gene == "nan" & !universe$gene == "" & !universe$gene == ".","gene"]


ggo1 <- enrichGO(gene = genes_clean, OrgDb = org.Hs.eg.db ,ont = "ALL", pvalueCutoff = 1, pAdjustMethod = "BH", universe = universe_clean, qvalueCutoff = 1, keyType = "SYMBOL",readable = F)

write.table(as.data.frame(ggo1), paste("cluster",element,"/go-enrichment-results-cluster",element,"-against-consensus.txt",sep=""))




}
 