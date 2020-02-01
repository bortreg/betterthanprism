#Load Packages
library(dplyr)
library(ggplot2)
library(tibble)

#Import Pathway Correlation Matrix
HUpaths <- read.csv("~/Desktop/HUpath.csv")
HEUpaths <- read.csv("~/Desktop/HEUpath.csv")
names(HUpaths) <- paste0("HU.", names(HUpaths))
names(HEUpaths) <- paste0("HEU.", names(HUpaths))
HUpathName <- HUpaths[,-1]
rownames(HUpathName) <- HUpaths[,1]
HEUpathName <- HEUpaths[,-1]
rownames(HEUpathName) <- HEUpaths[,1]

#Calculate Correlation Matrix
HUcormat <- signif(cor(HUpathName),2)
HEUcormat <- signif(cor(HEUpathName),2)

#Plot Correlation Matrix
col<- colorRampPalette(c("#4378EB", "#A09E9A", "#EEA826"))(20)
heatmap(HUcormat, col=col, symm=TRUE)
heatmap(HEUcormat, col=col, symm=TRUE)
