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

userUniverse <- readBed("all-consensus.bed")

### Making test set LOOP start

clustervector <- c("clusterI","clusterII","clusterIII","clusterIV","clusterV")


for (element in clustervector){
  
  print(element)
  
  userSets <- readBed(paste("",element,"-FC-clustered.bed",sep=""))
  
  ### Running LOLA
  
  locResults = runLOLA(userSets, userUniverse, regionDB, cores=1)
  
  writeCombinedEnrichment(locResults, outFolder = paste("all_consensus/",element,"/",sep=""), includeSplits=TRUE)
  
}

