#Load Packages
library(dplyr)
library(ggplot2)
library(tibble)
library(gplots)

#Import Pathway Correlation Matrix
HUpaths <- read.csv("~/Documents/R/Nanostring_DataFrames/HUpath.csv")
HEUpaths <- read.csv("~/Documents/R/Nanostring_DataFrames/HUpath.csv")
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
my_palette<- colorRampPalette(c("#FFFFFF", "#45A242"))(n = 40)
heatmap(HUcormat, col=my_palette, symm=TRUE)
heatmap(HEUcormat, col=my_palette, symm=TRUE)

##Heatmap of Normalized Reads##

#Load Normalized mRNA data
HUnorm <- read.csv("~/Documents/R/Nanostring_DataFrames/HUnorm.csv")
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
NFkBgenes <- c("TNF","IL1B","IL6","IL12A","IL12B","IL18","IL15","IFNG","IL4","IL10","IL1RAP","TGFB1","IL8","CCL2","CCL5","NFKB1","NFKB2","IL23A","condition")
HUnormNFKB <- mat_data[,NFkBgenes]
HUnormNFKB$cc <- ifelse(HUnormNFKB$condition == "BCG", "#ECC03F", "#C6C4BF")



#plot heatmap
pdf("~/Desktop/test.pdf", width=5, height=5)
heatmap <- heatmap.2(t(HUnormNFKB[,1:18]), 
           col=my_palette, 
           Rowv = FALSE, 
           dendrogram = "col",
           trace = "none",
           labCol = FALSE,
           margins = c(5,12),
           offsetRow = 0.01,
           keysize = 1.6,
           density.info = "none",
           key.par = list(cex=0.8),
           key.title = NA,
           key.xlab = "z-score",
           ColSideColors = HUnormNFKB$cc)

par(lend = 1)
legend("topright",
       legend = c("BCG", "untreated"),
       cex = 0.7,
       pt.cex = 1,
       col = c("#ECC03F", "#C6C4BF"),
       lty = 1,
       lwd = 10
)          
dev.off()

