library(ggplot2)
library(ggrepel)
library(plyr)
library(dplyr)
library(tidyverse)

# Importing Proteome Discoverer Protein Abundances. Accession numbers are converted to Row Names
pd.protein.abundance <- read.csv("CC_SKu_93-M790-B01-P9641-P9646-SMARCC1-Norm-20191127_Proteins.txt", head=T, sep="\t", dec = ".", na.strings = "NA", stringsAsFactors = F, row.names = "Accession") 

# Split Description column to multiple columns
pd.protein.abundance <- pd.protein.abundance %>% separate(Description, c("Description", "SV"), sep = "\\sSV=", remove = T, convert=T, extra = "merge", fill="right")

pd.protein.abundance <- pd.protein.abundance %>% separate(Description, c("Description", "PE"), sep = "\\sPE=", remove = T, convert=T, extra = "merge", fill="right")

pd.protein.abundance <- pd.protein.abundance %>% separate(Description, c("Description", "GN"), sep = "\\sGN=", remove = T, convert=T, extra = "merge", fill="right")

pd.protein.abundance <- pd.protein.abundance %>% separate(Description, c("Description", "OX"), sep = "\\sOX=", remove = T, convert=T, extra = "merge", fill="right")

pd.protein.abundance <- pd.protein.abundance %>% separate(Description, c("Description", "OS"), sep = "\\sOS=", remove = T, convert=T, extra = "merge", fill="right")

# Filter
pd.protein.abundance <- pd.protein.abundance[ which(pd.protein.abundance$Protein.FDR.Confidence.Combined=='High'), ]
pd.protein.abundance <- pd.protein.abundance[ which(pd.protein.abundance$X.Peptides > 1), ]

# Edit Column Names
colnames(pd.protein.abundance) <- str_replace_all(colnames(pd.protein.abundance), "\\.+", "\\.")

##########
# Output file
pdf(file="P9641_P9646.pdf", width=8.27, height=11.69)
# Output PDF , font and page dimensions
par(cex = 0.6);
par(mar = c(0,4,2,0))

# Comparison Set 1 (Chrom)
#############
wt.a4.cc1.ip.chrom <- pd.protein.abundance[, grepl("Abundances.Normalized.+Chrom", colnames(pd.protein.abundance))]
wt.a4.cc1.ip.chrom <- na.omit(wt.a4.cc1.ip.chrom)

wt.a4.cc1.ip.chrom <- merge(wt.a4.cc1.ip.chrom, pd.protein.abundance[c("GN", "Marked.as")], by=0, all.x=T)
#row.names(wt.a4.cc1.ip.chrom) <- wt.a4.cc1.ip.chrom$Row.names
#wt.a4.cc1.ip.chrom$Row.names <- NULL
wt.a4.cc1.ip.chrom$Marked.as <- gsub("BAF", T, wt.a4.cc1.ip.chrom$Marked.as)
wt.a4.cc1.ip.chrom$Marked.as <- as.logical(wt.a4.cc1.ip.chrom$Marked.as)
wt.a4.cc1.ip.chrom <- na.omit(wt.a4.cc1.ip.chrom)
wt.a4.cc1.ip.chrom$Marked.as <- NULL

wt.a4.cc1.ip.chrom[,2:5] <- log2(wt.a4.cc1.ip.chrom[,2:5])
colnames(wt.a4.cc1.ip.chrom) <- c("accession", paste("log2", colnames(wt.a4.cc1.ip.chrom[,2:5]), sep = "."), "gene.name")

wt.a4.cc1.ip.chrom$mean.dTAG.chrom <- apply(wt.a4.cc1.ip.chrom[,grepl("dTAG", colnames(wt.a4.cc1.ip.chrom))], 1, mean)
wt.a4.cc1.ip.chrom$mean.DMSO.chrom <- apply(wt.a4.cc1.ip.chrom[,grepl("DMSO", colnames(wt.a4.cc1.ip.chrom))], 1, mean)

wt.a4.cc1.ip.chrom$l2fc.chrom <- wt.a4.cc1.ip.chrom$mean.dTAG.chrom - wt.a4.cc1.ip.chrom$mean.DMSO.chrom

# Comparison Set 2 (Cyto)
#############
wt.a4.cc1.ip.cyto <- pd.protein.abundance[, grepl("Abundances.Normalized.+Cyto", colnames(pd.protein.abundance))]
wt.a4.cc1.ip.cyto <- na.omit(wt.a4.cc1.ip.cyto)

wt.a4.cc1.ip.cyto <- merge(wt.a4.cc1.ip.cyto, pd.protein.abundance[c("GN", "Marked.as")], by=0, all.x=T)
#row.names(wt.a4.cc1.ip.cyto) <- wt.a4.cc1.ip.cyto$Row.names
#wt.a4.cc1.ip.cyto$Row.names <- NULL
wt.a4.cc1.ip.cyto$Marked.as <- gsub("BAF", T, wt.a4.cc1.ip.cyto$Marked.as)
wt.a4.cc1.ip.cyto$Marked.as <- as.logical(wt.a4.cc1.ip.cyto$Marked.as)
wt.a4.cc1.ip.cyto <- na.omit(wt.a4.cc1.ip.cyto)
wt.a4.cc1.ip.cyto$Marked.as <- NULL

wt.a4.cc1.ip.cyto[,2:5] <- log2(wt.a4.cc1.ip.cyto[,2:5])
colnames(wt.a4.cc1.ip.cyto) <- c("accession", paste("log2", colnames(wt.a4.cc1.ip.cyto[,2:5]), sep = "."), "gene.name")

wt.a4.cc1.ip.cyto$mean.dTAG.cyto <- apply(wt.a4.cc1.ip.cyto[,grepl("dTAG", colnames(wt.a4.cc1.ip.cyto))], 1, mean)
wt.a4.cc1.ip.cyto$mean.DMSO.cyto <- apply(wt.a4.cc1.ip.cyto[,grepl("DMSO", colnames(wt.a4.cc1.ip.cyto))], 1, mean)

wt.a4.cc1.ip.cyto$l2fc.cyto <- wt.a4.cc1.ip.cyto$mean.dTAG.cyto - wt.a4.cc1.ip.cyto$mean.DMSO.cyto

# Comparison Set 3 (Nucl)
#############
wt.a4.cc1.ip.nucl <- pd.protein.abundance[, grepl("Abundances.Normalized.+Nucl", colnames(pd.protein.abundance))]
wt.a4.cc1.ip.nucl <- na.omit(wt.a4.cc1.ip.nucl)

wt.a4.cc1.ip.nucl <- merge(wt.a4.cc1.ip.nucl, pd.protein.abundance[c("GN", "Marked.as")], by=0, all.x=T)
#row.names(wt.a4.cc1.ip.nucl) <- wt.a4.cc1.ip.nucl$Row.names
#wt.a4.cc1.ip.nucl$Row.names <- NULL
wt.a4.cc1.ip.nucl$Marked.as <- gsub("BAF", T, wt.a4.cc1.ip.nucl$Marked.as)
wt.a4.cc1.ip.nucl$Marked.as <- as.logical(wt.a4.cc1.ip.nucl$Marked.as)
wt.a4.cc1.ip.nucl <- na.omit(wt.a4.cc1.ip.nucl)
wt.a4.cc1.ip.nucl$Marked.as <- NULL

wt.a4.cc1.ip.nucl[,2:5] <- log2(wt.a4.cc1.ip.nucl[,2:5])
colnames(wt.a4.cc1.ip.nucl) <- c("accession", paste("log2", colnames(wt.a4.cc1.ip.nucl[,2:5]), sep = "."), "gene.name")

wt.a4.cc1.ip.nucl$mean.dTAG.nucl <- apply(wt.a4.cc1.ip.nucl[,grepl("dTAG", colnames(wt.a4.cc1.ip.nucl))], 1, mean)
wt.a4.cc1.ip.nucl$mean.DMSO.nucl <- apply(wt.a4.cc1.ip.nucl[,grepl("DMSO", colnames(wt.a4.cc1.ip.nucl))], 1, mean)

wt.a4.cc1.ip.nucl$l2fc.nucl <- wt.a4.cc1.ip.nucl$mean.dTAG.nucl - wt.a4.cc1.ip.nucl$mean.DMSO.nucl

####
library(pheatmap)
library(grid)

wt.a4.cc1.ip.combined <- merge(wt.a4.cc1.ip.cyto, wt.a4.cc1.ip.nucl, by=intersect(colnames(wt.a4.cc1.ip.cyto), colnames(wt.a4.cc1.ip.nucl)), all=T, sort=F)
wt.a4.cc1.ip.combined <- merge(wt.a4.cc1.ip.combined, wt.a4.cc1.ip.chrom, by=intersect(colnames(wt.a4.cc1.ip.combined), colnames(wt.a4.cc1.ip.chrom)), all=T, sort=F)
row.names(wt.a4.cc1.ip.combined) <- wt.a4.cc1.ip.combined$gene.name
wt.a4.cc1.ip.combined[is.na(wt.a4.cc1.ip.combined)] <- 0

wt.subset.mean <- wt.a4.cc1.ip.combined[c("mean.DMSO.cyto", "mean.dTAG.cyto", "mean.DMSO.nucl", "mean.dTAG.nucl", "mean.DMSO.chrom", "mean.dTAG.chrom")]
# wt.a4.cc1.ip.combined[,grepl("mean", colnames(df1))]

pheatmap(wt.subset.mean, main= "mean.log2.Normalized.Abundances", cluster_cols=F, fontsize_row = 20, show_rownames = T, fontsize=8, labels_col=c("Cyto DMSO","Cyto dTAG","Nucl DMSO","Nucl dTAG","Chrom DMSO","Chrom dTAG"), fontsize_number=10, fontsize_col = 10, cellwidth = 60, color = colorRampPalette(c("white","yellow", "blue"))(50))

wt.subset.l2fc <- wt.a4.cc1.ip.combined[c("l2fc.cyto", "l2fc.nucl", "l2fc.chrom")]
wt.subset.l2fc.scale <- max ( abs(min(wt.subset.l2fc)), max(wt.subset.l2fc))

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
#myBreaks <- c(seq(min(df1[,grepl("l2fc", colnames(df1))]), 0, length.out=ceiling(paletteLength/2) + 1), 
#              seq(max(df1[,grepl("l2fc", colnames(df1))])/paletteLength, max(df1[,grepl("l2fc", colnames(df1))]), length.out=floor(paletteLength/2)))

myBreaks <- c(seq(-(wt.subset.l2fc.scale), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(wt.subset.l2fc.scale/paletteLength, wt.subset.l2fc.scale, length.out=floor(paletteLength/2)))

pheatmap(wt.subset.l2fc, main= "log2.fold.change.Normalized.Abundances", cluster_cols=F, fontsize_row = 20, show_rownames = T, fontsize=8, labels_col=c("Cyto","Nucl", "Chrom"), fontsize_number=10, fontsize_col = 10, cellwidth = 60, color=myColor, breaks=myBreaks)

dev.off()

row.names(wt.a4.cc1.ip.combined) <- wt.a4.cc1.ip.combined$accession
wt.a4.cc1.ip.pd.combined <- merge(wt.a4.cc1.ip.combined, pd.protein.abundance, by=0, all.x=T)
wt.a4.cc1.ip.pd.combined$Row.names <- NULL
write.table(wt.a4.cc1.ip.pd.combined, file="P9641_P9646.csv", row.names = F, sep="\t", quote = F)

#pd.protein.abundance2 <- merge(pd.protein.ab undance, wt.a4.cc1.ip.chrom, by=0, all.x=T)
#row.names(pd.protein.abundance2) <- pd.protein.abundance2$Row.names
#pd.protein.abundance2$Row.names <- NULL
