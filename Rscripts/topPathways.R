library(dplyr)
library(gplots)

#Import Pathway Scores
pathScores <- read.csv("~/Desktop/pathMerge.csv")
pathScores <- pathScores[order(pathScores$HU.MeanBCG, decreasing = TRUE),]

Pathtop50 <- pathScores[1:50, c("HU.MeanUntreat","HU.MeanBCG")] 
row.names(Pathtop50) <- pathScores[1:50,"Pathway"]
Pathtop50 <- t(Pathtop50)  


#Heatmap of path score top 50
my_palette<- colorRampPalette(c("#E2ECE1","#C6E0C5","#89BE88","#5CA25A","#45A242"))(n = 40.0)
pdf("~/Desktop/test.pdf", width=5, height = 10)
heatmap <- heatmap.2(t(Pathtop50[,1:50]), 
                     col=my_palette, 
                     Rowv = FALSE, 
                     dendrogram = "none",
                     trace = "none",
                     labCol = FALSE,
                     margins = c(10,25),
                     density.info = "none",
                     Colv = FALSE,
                     offsetRow = 0.01,
                     key.par = list(cex=0.8),
                     key.title = NA,
                     key.xlab = "pathway scores",
                     ColSideColors = c("#ECC03F", "#C6C4BF"))

par(lend = 1)
legend("topright",
       legend = c("Untreated", "BCG"),
       cex = 0.7,
       pt.cex = 1,
       col = c("#ECC03F", "#C6C4BF"),
       lty = 1,
       lwd = 10
)          
dev.off()
