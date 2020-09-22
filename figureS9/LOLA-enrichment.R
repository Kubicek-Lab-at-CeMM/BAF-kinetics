library(LOLA)
library(simpleCache)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(RColorBrewer)

setwd("../data/")

### setting database

dbPath = "LOLA_db/hg_38/"
regionDB = loadRegionDB(dbPath)

### analyses consensus

### Reading in universe set

userUniverse <- readBed("universe-consensus.bed")

### Making test set LOOP start

clustervector <- c("cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8","cluster9","cluster10","cluster11")


for (element in clustervector){
  
  print(element)
  
  userSets <- readBed(paste("",element,".bed",sep=""))
  
  ### Running LOLA
  
  locResults = runLOLA(userSets, userUniverse, regionDB, cores=1)
  
  writeCombinedEnrichment(locResults, outFolder = paste("all_consensus/",element,"/",sep=""), includeSplits=TRUE)
  
}

