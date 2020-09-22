### Getting packages

library(data.table)
library(DESeq2)

### setting directory

setwd("../data/")

### reading data 

### WT SMARCA4 dtag dtag 3h

  
  A <- read.table(paste("WT_SMARCA4_dtag_3h_dtag_A_Ensembl_genebody_unnorm.cov",sep = ""),header = F)
  
  colnames(A) <- c("chr","pos1","pos2","gene","ensembl_id","strand",paste("3h_dtag_WT_SMARCA4_A",sep = ""))
  
  B <- read.table(paste("WT_SMARCA4_dtag_3h_dtag_B_Ensembl_genebody_unnorm.cov",sep = ""),header = F)
  
  colnames(B) <- c("chr","pos1","pos2","gene","ensembl_id","strand",paste("3h_dtag_WT_SMARCA4_B",sep = ""))
  
  A_DMSO <- read.table("WT_SMARCA4_dtag_3h_DMSO_A_Ensembl_genebody_unnorm.cov",header = F)
  
  colnames(A_DMSO) <- c("chr","pos1","pos2","gene","ensembl_id","strand","3h_DMSO_WT_SMARCA4_A")
  
  B_DMSO <- read.table("WT_SMARCA4_dtag_3h_DMSO_B_Ensembl_genebody_unnorm.cov",header = F)
  
  colnames(B_DMSO) <- c("chr","pos1","pos2","gene","ensembl_id","strand","3h_DMSO_WT_SMARCA4_B")
  
  countdata <- data.frame(A[,7],B[,7],A_DMSO[,7],B_DMSO[,7])
  
  colnames(countdata) <- c(paste("3h_dtag_WT_SMARCA4_A",sep = ""),paste("3h_dtag_WT_SMARCA4_B",sep = ""),"3h_DMSO_WT_SMARCA4_A","3h_DMSO_WT_SMARCA4_B")
  
  row.names(countdata) <- A$ensembl_id
  
  coldata <- data.frame(t(data.frame("3h_dtag","3h_dtag","DMSO","DMSO")))
  
  row.names(coldata) <- colnames(countdata)
  
  colnames(coldata) <- "condition"
  
  ### Generate DESeq object 
  
  dds <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = coldata,
                                design = ~ condition)
  
  
  dds$condition <- relevel(dds$condition, ref = "DMSO")
  
  analysis_res <- DESeq(dds)
  
  results_res <- results(analysis_res)
  
  results_res$comparison <- paste("3h_dtag_WT_SMARCA4_dtag_vs_3h_DMSO_WT_SMARCA4_dtag",sep = "")
  
  results_res$index <- paste(A$chr,":",A$pos1,"-",A$pos2,sep="")
  
  results_res$gene <- A$gene
  
  results_res$ensembl_id <- A$ensembl_id
  
  
  write.table(results_res, paste("deseq2-3h-dtag-WT-SMARCA4-dtag-vs-3h-DMSO-WT-SMARCA4-dtag-gene-full-genebody.txt",sep = ""),quote=F,sep="\t",row.names=F)
  
  