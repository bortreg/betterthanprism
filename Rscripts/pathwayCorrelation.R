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
my_palette<- colorRampPalette(c("#FFFFFF", "#E89D11"))(n = 299)
heatmap(HUcormat, col=my_palette, symm=TRUE)
heatmap(HEUcormat, col=my_palette, symm=TRUE)

#Load Normalized mRNA data
HUnorm <- read.csv("~/Desktop/HUnorm.csv")
HUnorm <- as_tibble(HUnorm)

#Prepare data for heatmap
zscore <- function(x) {          #zscore function
  z <- (mat_data - mean(x)) / sd(x)
  return(z)
}
rnames <- as.matrix(HUnorm[,1])
mat_data <- data.matrix(HUnorm[,4:ncol(HUnorm)])
rownames(mat_data) <- rnames
mat_data <- log2(mat_data)
mat_data <- zscore(mat_data)
mat_data <- cbind(mat_data, HUnorm[,3])

#Index Genes
NFkBgenes <- c("TNF.mRNA","IL1B.mRNA","IL6.mRNA","IL12A.mRNA","IL12B.mRNA","IL18.mRNA","IL15.mRNA","IFNG.mRNA","IL4.mRNA","IL10.mRNA","IL1RAP.mRNA","TGFB1.mRNA","IL8.mRNA","CCL2.mRNA","CCL5.mRNA")
HUnormNFKB <- mat_data[,NFkBgenes]

#plot heatmap
heatmap.2(t(HUnormNFKB), 
          col=my_palette, 
          Rowv = FALSE, 
          dendrogram = "col",
          trace = "none",
          ColSideColors = c(
            rep("")
          

