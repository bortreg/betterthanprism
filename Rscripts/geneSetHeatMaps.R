library(dplyr)
library(RColorBrewer)
library(gplots)

##Produce Heat Maps for Each Gene Set in 2 Groups

#set color palette
my_palette<- colorRampPalette(c("black","gray","red"))(n = 40)

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

#Load Gene Sets
geneSets <- read.csv("~/Documents/R/Nanostring_DataFrames/GeneModuleLists.csv")
geneSets <- as_tibble(geneSets)

#Define
inflammGenes <- c(as.character(geneSets$inflammatory.response[1:35]), "condition")
HUnormInflamm <- mat_data[,inflammGenes]
HUnormInflamm$cc <- ifelse(HUnormInflamm$condition == "BCG", "#ECC03F", "#2E2EFE")

apopGenes <- c()
HUnormApop <- mat_data[,apopGenes]
HUnormApop$cc <- ifelse(HUnormApop$condition == "BCG", "#ECC03F", "#2E2EFE")

jakStatGenes <- c()
HUnormJakStat <- mat_data[,jakStatGenes]
HUnormJakStat$cc <- ifelse(HUnormJakStat$condition == "BCG", "#ECC03F", "#2E2EFE")

cellProlifGenes <- c()
HUnormCellProlif <- mat_data[,cellProlifGenes]
HUnormCellProlif$cc <- ifelse(HUnormCellProlif$condition == "BCG", "#ECC03F", "#2E2EFE")

nkbCasGenes <- c()
HUnormNkbCas <- mat_data[,nkbCasGenes]
HUnormNkbCas$cc <- ifelse(HUnormNkbCas$condition == "BCG", "#ECC03F", "#2E2EFE")

#plot heatmap
hm_df <- HUnormInflamm
pdf("~/Desktop/test2.pdf", width=5, height=((length(hm_df)-2)*0.285))
heatmap <- heatmap.2(t(hm_df[,1:(length(hm_df)-2)]), 
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
                     ColSideColors = HUnormMycoB$cc)

par(lend = 1)
legend("topright",
       legend = c("untreated","BCG"),
       cex = 0.7,
       pt.cex = 1,
       col = c("#2E2EFE","#ECC03F"),
       lty = 1,
       lwd = 10
)          
dev.off()
