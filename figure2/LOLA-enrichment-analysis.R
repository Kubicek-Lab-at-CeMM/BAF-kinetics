library(LOLA)
library(simpleCache)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(RColorBrewer)

setwd("/Users/sgrosche/Development/191127-BAF-timecourse-batch3/results/differential_analysis/191206-SMARCA4-timecourse-and-KO-FC-based-3-clusters/200129_LOLA_enrichment/")

### setting database

dbPath = "/Users/sgrosche/Development/190917-resources/LOLA_db/hg_38/"
regionDB = loadRegionDB(dbPath)

### analyses consensus

### Reading in universe set

userUniverse <- readBed("/Users/sgrosche/Development/191127-BAF-timecourse-batch3/results/differential_analysis/191206-SMARCA4-timecourse-and-KO-FC-based-3-clusters/LOLA_enrichment/191209-universe-consensus.bed")

### Making test set LOOP start

clustervector <- c("cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8","cluster9","cluster10")


for (element in clustervector){
  
  print(element)
  
  userSets <- readBed(paste("/Users/sgrosche/Development/200128-cluster-regions/",element,".bed",sep=""))
  
  ### Running LOLA
  
  locResults = runLOLA(userSets, userUniverse, regionDB, cores=1)
  
  writeCombinedEnrichment(locResults, outFolder = paste("all_consensus/",element,"/",sep=""), includeSplits=TRUE)
  
}