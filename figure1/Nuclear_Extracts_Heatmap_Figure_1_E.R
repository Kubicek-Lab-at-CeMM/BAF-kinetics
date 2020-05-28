LibsRequired = c("plyr", "dplyr", "tidyverse", "ggplot2", "ggrepel", "stringr", "RColorBrewer")

# Select - WTA4, CC1KOCC2, A4KOCA2, A2KOA4
# Importing Proteome Discoverer Protein Abundances. Accession numbers are converted to Row Names
pd.protein.abundance <- read.csv("~/Documents/Sandra/CC_SKu_93-M779-B07-P9510_Proteins.txt", head=T, sep="\t", dec = ".", na.strings = "NA", stringsAsFactors = F, row.names = "Accession") 

pd.protein.abundance <- read.csv("/home/nmarella/Documents/Sandra/CC_SKu_93-CC1KOCC2-201911_Proteins.txt", head=T, sep="\t", dec = ".", na.strings = "NA", stringsAsFactors = F, row.names = "Accession") 

pd.protein.abundance <- read.csv("/home/nmarella/Documents/Sandra/CC_SKu_93-A4KOCA2-201911_Proteins.txt", head=T, sep="\t", dec = ".", na.strings = "NA", stringsAsFactors = F, row.names = "Accession") 

pd.protein.abundance <- read.csv("/home/nmarella/Documents/Sandra/CC_SKu_93-A2KOA4-20191125_Proteins.txt", head=T, sep="\t", dec = ".", na.strings = "NA", stringsAsFactors = F, row.names = "Accession") 

# For all- WTA4, CC1KOCC2, A4KOCA2, A2KOA4
# Split Description column to multiple columns
pd.protein.abundance <- pd.protein.abundance %>% separate(Description, c("Description", "SV"), sep = "\\sSV=", remove = T, convert=T, extra = "merge", fill="right")

# For all- WTA4, CC1KOCC2, A4KOCA2, A2KOA4
pd.protein.abundance <- pd.protein.abundance %>% separate(Description, c("Description", "PE"), sep = "\\sPE=", remove = T, convert=T, extra = "merge", fill="right")

# For all- WTA4, CC1KOCC2, A4KOCA2, A2KOA4
pd.protein.abundance <- pd.protein.abundance %>% separate(Description, c("Description", "GN"), sep = "\\sGN=", remove = T, convert=T, extra = "merge", fill="right")

# For all- WTA4, CC1KOCC2, A4KOCA2, A2KOA4
pd.protein.abundance <- pd.protein.abundance %>% separate(Description, c("Description", "OX"), sep = "\\sOX=", remove = T, convert=T, extra = "merge", fill="right")

# For all- WTA4, CC1KOCC2, A4KOCA2, A2KOA4
pd.protein.abundance <- pd.protein.abundance %>% separate(Description, c("Description", "OS"), sep = "\\sOS=", remove = T, convert=T, extra = "merge", fill="right")

# Filter
pd.protein.abundance <- pd.protein.abundance[ which(pd.protein.abundance$Protein.FDR.Confidence.Combined=='High'), ]
pd.protein.abundance <- pd.protein.abundance[ which(pd.protein.abundance$Number.of.Peptides > 1), ]

# Edit Column Names
# For all- WTA4, CC1KOCC2, A4KOCA2, A2KOA4
colnames(pd.protein.abundance) <- str_replace_all(colnames(pd.protein.abundance), "\\.+", "\\.")

# T Test Function
myTtest <- function(df, grp1, grp2) {
     x = df[grp1]
     y = df[grp2]
     x = as.numeric(x)
     y = as.numeric(y)  
     results = t.test(x, y)
     results$p.value
}

my.sum <- function(x) {
  sum(x[!is.na(x)])
}

# Select For WTA4
# only BAF  
baf.pd <- pd.protein.abundance[ which(pd.protein.abundance$Marked.as=='BAF'),]
baf.pd <- baf.pd[, grepl("Abundances.Normalized", colnames(baf.pd))]

# For WTA4
baf.pd$dmso.bio.rep.1 <- apply( baf.pd[, grepl("DMSO\\.1\\.1|DMSO\\.1\\.2", colnames(baf.pd))] , 1, my.sum )
baf.pd$dmso.bio.rep.2 <- apply( baf.pd[, grepl("DMSO\\.2\\.1|DMSO\\.2\\.2", colnames(baf.pd))] , 1, my.sum )

# For WTA4
baf.pd$dtag.72h.bio.rep.1 <- apply( baf.pd[, grepl("dTAG\\.72h\\.1\\.1|dTAG\\.72h\\.1\\.2", colnames(baf.pd))] , 1, my.sum )
baf.pd$dtag.72h.bio.rep.2 <- apply( baf.pd[, grepl("dTAG\\.72h\\.2\\.1|dTAG\\.72h\\.2\\.2", colnames(baf.pd))] , 1, my.sum )


# For WTA4
baf.pd$dtag.24h.bio.rep.1 <- apply( baf.pd[, grepl("dTAG\\.24h\\.1\\.3|dTAG\\.24h\\.1\\.4", colnames(baf.pd))] , 1, my.sum )
baf.pd$dtag.24h.bio.rep.2 <- apply( baf.pd[, grepl("dTAG\\.24h\\.2\\.3|dTAG\\.24h\\.2\\.4", colnames(baf.pd))] , 1, my.sum )

# For WTA4
baf.pd$dtag.6h.bio.rep.1 <- apply( baf.pd[, grepl("dTAG\\.6h\\.1\\.1|dTAG\\.6h\\.1\\.2", colnames(baf.pd))] , 1, my.sum )
baf.pd$dtag.6h.bio.rep.2 <- apply( baf.pd[, grepl("dTAG\\.6h\\.2\\.1|dTAG\\.6h\\.2\\.2", colnames(baf.pd))] , 1, my.sum )

# For WTA4
baf.pd$dtag.3h.bio.rep.1 <- apply( baf.pd[, grepl("dTAG\\.3h\\.1\\.1|dTAG\\.3h\\.1\\.2", colnames(baf.pd))] , 1, my.sum )
baf.pd$dtag.3h.bio.rep.2 <- apply( baf.pd[, grepl("dTAG\\.3h\\.2\\.1|dTAG\\.3h\\.2\\.2", colnames(baf.pd))] , 1, my.sum )

# For WTA4
baf.pd$dtag.2h.bio.rep.1 <- apply( baf.pd[, grepl("dTAG\\.2h\\.1\\.1|dTAG\\.2h\\.1\\.2", colnames(baf.pd))] , 1, my.sum )
baf.pd$dtag.2h.bio.rep.2 <- apply( baf.pd[, grepl("dTAG\\.2h\\.2\\.1|dTAG\\.2h\\.2\\.2", colnames(baf.pd))] , 1, my.sum )

# For WTA4, CC1KOCC2, A4KOCA2, A2KOA4
baf.pd <-  baf.pd[, !grepl("Abundances.Normalized", colnames(baf.pd))]
#baf.pd <- na.omit(baf.pd) 
baf.pd <- log2(baf.pd)
colnames(baf.pd) <- c(paste("log2", colnames(baf.pd), sep = "."))
is.na(baf.pd) <- do.call(cbind,lapply(baf.pd, is.infinite))
baf.pd <- na.omit(baf.pd)
baf.pd$log2.mean.dmso <-  apply( baf.pd[, grepl("dmso", colnames(baf.pd))] , 1, mean )
baf.pd$log2.mean.72 <-  apply( baf.pd[, grepl("dtag\\.72h", colnames(baf.pd))] , 1, mean )

# For WTA4
baf.pd$log2.mean.24 <-  apply( baf.pd[, grepl("dtag\\.24h", colnames(baf.pd))] , 1, mean )
baf.pd$log2.mean.6 <-  apply( baf.pd[, grepl("dtag\\.6h", colnames(baf.pd))] , 1, mean )
baf.pd$log2.mean.3 <-  apply( baf.pd[, grepl("dtag\\.3h", colnames(baf.pd))] , 1, mean )
baf.pd$log2.mean.2 <-  apply( baf.pd[, grepl("dtag\\.2h", colnames(baf.pd))] , 1, mean )

# For WTA4, CC1KOCC2, A4KOCA2, A2KOA4
baf.pd$l2fc.72h <- baf.pd$log2.mean.72 - baf.pd$log2.mean.dmso 
baf.pd$p.72h = apply(baf.pd, 1, myTtest, grp1 = grep("dtag\\.72h\\.bio", colnames(baf.pd), value=T), grp2 = grep("dmso\\.bio", colnames(baf.pd), value=T))
baf.pd$p.adj.72h <- p.adjust(baf.pd$p.72h, method = "fdr")

# For WTA4
baf.pd$l2fc.24h <- baf.pd$log2.mean.24 - baf.pd$log2.mean.dmso 
baf.pd$p.24h = apply(baf.pd, 1, myTtest, grp1 = grep("dtag\\.24h\\.bio", colnames(baf.pd), value=T), grp2 = grep("dmso\\.bio", colnames(baf.pd), value=T))
baf.pd$p.adj.24h <- p.adjust(baf.pd$p.24h, method = "fdr")

# For WTA4
baf.pd$l2fc.6h <- baf.pd$log2.mean.6 - baf.pd$log2.mean.dmso 
baf.pd$p.6h = apply(baf.pd, 1, myTtest, grp1 = grep("dtag\\.6h\\.bio", colnames(baf.pd), value=T), grp2 = grep("dmso\\.bio", colnames(baf.pd), value=T))
baf.pd$p.adj.6h <- p.adjust(baf.pd$p.6h, method = "fdr")

# For WTA4
baf.pd$l2fc.3h <- baf.pd$log2.mean.3 - baf.pd$log2.mean.dmso 
baf.pd$p.3h = apply(baf.pd, 1, myTtest, grp1 = grep("dtag\\.3h\\.bio", colnames(baf.pd), value=T), grp2 = grep("dmso\\.bio", colnames(baf.pd), value=T))
baf.pd$p.adj.3h <- p.adjust(baf.pd$p.3h, method = "fdr")

# For WTA4
baf.pd$l2fc.2h <- baf.pd$log2.mean.2 - baf.pd$log2.mean.dmso 
baf.pd$p.2h = apply(baf.pd, 1, myTtest, grp1 = grep("dtag\\.2h\\.bio", colnames(baf.pd), value=T), grp2 = grep("dmso\\.bio", colnames(baf.pd), value=T))
baf.pd$p.adj.2h <- p.adjust(baf.pd$p.2h, method = "fdr")

# For WTA4, CC1KOCC2, A4KOCA2, A2KOA4
baf.pd <- merge(baf.pd, pd.protein.abundance, by=0, all.x=T)
baf.pd$GN[which(is.na(baf.pd$GN))] <- baf.pd$Row.names[which(is.na(baf.pd$GN))]
baf.pd$GN[which(duplicated(baf.pd$GN))] <- baf.pd$Row.names[which(duplicated(baf.pd$GN))]
row.names(baf.pd) <- baf.pd$GN

# Select - WTA4, CC1KOCC2, A4KOCA2, A2KOA4
write.table(baf.pd, file="~/Documents/Sandra/P9510_P9537_Nuclear_WTA4_TimeCourse.csv", row.names = F, sep="\t", quote = F)
write.table(baf.pd, file="~/Documents/Sandra/P9510_P9537_Nuclear_CC1KOCC2.csv", row.names = F, sep="\t", quote = F)
write.table(baf.pd, file="~/Documents/Sandra/P9510_P9537_Nuclear_A4KOCA2.csv", row.names = F, sep="\t", quote = F)
write.table(baf.pd, file="~/Documents/Sandra/P9510_P9537_Nuclear_A2KOA4.csv", row.names = F, sep="\t", quote = F)
      
  # For WTA4
  # Heatmap log 2 fold change
  
      baf.pd$l2fc.72h[which(baf.pd$p.72h > 0.05)] <- 0
      baf.pd$l2fc.24h[which(baf.pd$p.24h > 0.05)] <- 0
      baf.pd$l2fc.6h[which(baf.pd$p.6h > 0.05)] <- 0
      baf.pd$l2fc.3h[which(baf.pd$p.3h > 0.05)] <- 0
      baf.pd$l2fc.2h[which(baf.pd$p.2h > 0.05)] <- 0
  
  baf.pd.l2fc <-  baf.pd[c("l2fc.2h","l2fc.3h","l2fc.6h","l2fc.24h","l2fc.72h")]
  paletteLength <- 50
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  baf.pd.l2fc.scale <- max ( abs(min(baf.pd.l2fc)), max(baf.pd.l2fc))
  myBreaks <- c(seq(-(baf.pd.l2fc.scale), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(baf.pd.l2fc.scale/paletteLength, baf.pd.l2fc.scale, length.out=floor(paletteLength/2)))
  #myBreaks <- c(seq(min(baf.pd.l2fc), 0, length.out=ceiling(paletteLength/2) + 1), seq(max(baf.pd.l2fc)/paletteLength, max(baf.pd.l2fc), length.out=floor(paletteLength/2)))

  # Only BAF  
  pdf(file="~/Documents/Sandra/P9510_P9537_Nuclear_WTA4_TimeCourse_BAF_Heatmap.pdf", width=8.27, height=11.69)
  par(cex = 0.6)
  par(mar = c(0,4,2,0))
  pheatmap(baf.pd.l2fc, main= "BAF log2 fold change", cluster_cols=F, fontsize_row = 10, show_rownames = T, fontsize=10, labels_col=c("2h","3h", "6h", "24h", "72h"), fontsize_number=10, fontsize_col = 20, cellwidth = 60, color=myColor, breaks=myBreaks, angle_col=45, display_numbers = F)
  
  dev.off()
