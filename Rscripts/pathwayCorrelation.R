#Load Packages
library(dplyr)
library(ggplot2)
library(tibble)
library(gplots)

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
my_palette<- colorRampPalette(c("#4378EB", "#A09E9A", "#EEA826"))(20)
heatmap(HUcormat, col=my_palette, symm=TRUE)
heatmap(HEUcormat, col=my_palette, symm=TRUE)

#Load Normalized mRNA data
HUnorm <- read.csv("~/Desktop/HUnorm.csv")
HUnormName <- HUnorm[,-1]
rownames(HUnormName) <- HUnorm[,1]

#Transform Data
HUnormName <- log2(HUnormName)

HUnormMat <- as.matrix(HUnormName)
heatmap(HUnormMat, col=my_palette, scale = "none")
