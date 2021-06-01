library(dplyr)
library(RColorBrewer)
library(gplots)
library(ggplot2)

##Produce Heat Maps for Each Gene Set in 2 Groups

#set color palette
my_palette<- colorRampPalette(c("white","gray","#57C84D"))(n = 40)

#Load Normalized mRNA data
HInorm <- read.csv("~/Documents/R/Nanostring_DataFrames/HInorm.csv")
HInorm <- as_tibble(HInorm)
HEUnorm <- read.csv("~/Documents/R/Nanostring_DataFrames/HEUnorm.csv")
HEUnorm <- as_tibble(HEUnorm)

#Prepare data for heatmap, otherwise skip to barplots 
zscore <- function(x) {          #zscore function
  z <- (mat_data - mean(x)) / sd(x)
  return(z)
}
rnames <- as.matrix(HInorm[,1])
mat_data <- data.matrix(HInorm[,4:ncol(HInorm)])
rownames(mat_data) <- rnames
mat_data <- log2(mat_data)
mat_data <- zscore(mat_data)
mat_data <- cbind(mat_data, HInorm[,3])

#Load Condensed Gene Set List
geneSets <- read.csv("~/Documents/R/Nanostring_DataFrames/GeneModuleCond2.csv")
geneSets <- as_tibble(geneSets)
condGenes <- c(as.character(geneSets$condensed), "condition")
HInormCond <- mat_data
HInormCond$cc <- ifelse(HInormCond$condition == "SARS2uv", "#ECC03F", "#2E2EFE")

#plot heatmap
hm_df <- HInormCond[,1:50]
pdf("~/Desktop/test3.pdf", width=8, height=((length(hm_df)-2)*0.18))
heatmap <- heatmap.2(t(hm_df[,1:(length(hm_df)-2)]), 
                     col=my_palette, 
                     Rowv = TRUE, 
                     dendrogram = "row",
                     trace = "none",
                     labCol = FALSE,
                     margins = c(7,12),
                     offsetRow = 0.01,
                     keysize = 1.8, 
                     density.info = "none",
                     key.par = list(cex=0.8),
                     key.title = NA,
                     key.xlab = "z-score",
                     ColSideColors = HInormCond$cc,
                     lhei = c(1,10))
par(lend = 1)
legend("topright",
       legend = c("untreated","SARS2uv"),
       cex = 0.7,
       pt.cex = 1,
       col = c("#2E2EFE","#ECC03F"),
       lty = 1,
       lwd = 10
)             
dev.off()

#Use barplots to compare norm mRNA counts for individual genes
HUrnames <- as.matrix(HUnorm[,1])
HUcounts <- data.matrix(HUnorm[,4:ncol(HUnorm)])
rownames(HUcounts) <- HUrnames
HUcounts <- cbind(HUcounts, HUnorm[,3])
HUcounts <- HUcounts[order(HUcounts$condition, decreasing = FALSE),]
ggplot(HUcounts, aes(HUcounts$condition, HUcounts$CD40)) +
geom_bar(stat = "identity") +
labs(title = NULL, x = "Condition", y = "Normalized mRNA counts")

HEUrnames <- as.matrix(HEUnorm[,1])
HEUcounts <- data.matrix(HEUnorm[,4:ncol(HEUnorm)])
rownames(HEUcounts) <- HEUrnames
HEUcounts <- cbind(HEUcounts, HEUnorm[,3])
HEUcounts <- HEUcounts[order(HEUcounts$condition, decreasing = FALSE),]
ggplot(HEUcounts, aes(HEUcounts$condition, HEUcounts$CD40)) +
  geom_bar(stat = "identity") +
  labs(title = NULL, x = "Condition", y = "Normalized mRNA counts")



#Define Gene Sets
inflammGenes <- c(as.character(geneSets$inflammatory.response[1:35]), "condition")
HUnormInflamm <- mat_data[,inflammGenes]
HUnormInflamm$cc <- ifelse(HUnormInflamm$condition == "BCG", "#ECC03F", "#2E2EFE")

apopGenes <- c(as.character(geneSets$programmed.cell.death[1:49]), "condition")
HUnormApop <- mat_data[,apopGenes]
HUnormApop$cc <- ifelse(HUnormApop$condition == "BCG", "#ECC03F", "#2E2EFE")

jakStatGenes <- c(as.character(geneSets$jak.stat.cascade[1:15]), "condition")
HUnormJakStat <- mat_data[,jakStatGenes]
HUnormJakStat$cc <- ifelse(HUnormJakStat$condition == "BCG", "#ECC03F", "#2E2EFE")

cellProlifGenes <- c(as.character(geneSets$cell.proliferation[1:19]), "condition")
HUnormCellProlif <- mat_data[,cellProlifGenes]
HUnormCellProlif$cc <- ifelse(HUnormCellProlif$condition == "BCG", "#ECC03F", "#2E2EFE")

nkbCasGenes <- c(as.character(geneSets$Ikk.nfkb.cascade[1:20]), "condition")
HUnormNkbCas <- mat_data[,nkbCasGenes]
HUnormNkbCas$cc <- ifelse(HUnormNkbCas$condition == "BCG", "#ECC03F", "#2E2EFE")
  