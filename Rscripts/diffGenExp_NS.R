#Load Packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(RColorBrewer)
library(gplots)

#Import Tables of Differentially Expressed Genes from nSolver
#*Make sure to rename columns you wish to remain unique to a specific group
HU_DE <- read.csv("~/Documents/R/Nanostring_DataFrames/DE_HU.csv")
HEU_DE <- read.csv("~/Documents/R/Nanostring_DataFrames/DE_HEU.csv")

#Plot volcanoes
df <- HU_DE
pdf("~/Desktop/test.pdf", width=5, height=5)
ggplot(df, aes(x=Log2.fold.change, y=-log2(P.value), color = FDR)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("#45A242","#5CA25A","#89BE88","#C6E0C5","#CECFCE")) +
  theme_bw() +
  geom_text_repel(
    data = df[1:50,],
    aes(label = Gene),
    color = "black",
    size = 2,
    segment.size = 0.1,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  geom_hline(yintercept = -log2(0.000181), linetype="solid", color = "darkgray", size=0.5) +
  geom_hline(yintercept = -log2(0.00129), linetype="dotted", color = "darkgray", size=0.5) +
  geom_hline(yintercept = -log2(0.00294), linetype="dashed", color = "darkgray", size=0.5) +
  geom_hline(yintercept = -log2(0.0251), linetype="dotdash", color = "darkgray", size=0.5)
dev.off()

#Merge dataframes
names(HU_DE) <- paste0("HU.", names(HU_DE))
names(HU_DE) <- gsub("HU.Gene","Gene",names(HU_DE))
names(HEU_DE) <- paste0("HEU.", names(HEU_DE))
names(HEU_DE) <- gsub("HEU.Gene","Gene",names(HEU_DE))

DEcomb <- merge(HU_DE, HEU_DE)
DEcomb <- DEcomb[order(DEcomb$HU.Log2.fold.change, decreasing = TRUE), ]
DEtop50 <- DEcomb[1:50, c("HU.Linear.fold.change","HEU.Linear.fold.change")] 
row.names(DEtop50) <- DEcomb[1:50,"Gene"]
DEtop50 <- t(DEtop50)  

#Label DE gene hits



#Heatmap of DEtop100
DEtop50 <- log2(DEtop50)
my_palette<- colorRampPalette(c("#E2ECE1","#C6E0C5","#89BE88","#5CA25A","#45A242"))(n = 40.0)
pdf("~/Desktop/test.pdf", width=3, height = 8)
heatmap <- heatmap.2(t(DEtop50[,1:50]), 
                     col=my_palette, 
                     Rowv = FALSE, 
                     dendrogram = "none",
                     trace = "none",
                     labCol = FALSE,
                     margins = c(5,5),
                     density.info = "none",
                     Colv = FALSE,
                     offsetRow = 0.01,
                     key.par = list(cex=0.8),
                     key.title = NA,
                     key.xlab = "log2 Fold Change",
                     ColSideColors = c("#838284", "#9554BE"))

par(lend = 1)
legend("topright",
       legend = c("HU", "HEU"),
       cex = 0.7,
       pt.cex = 1,
       col = c("#838284", "#9554BE"),
       lty = 1,
       lwd = 10
)          
dev.off()


#Calculate z scores from ratios of fold change from experimental group over control group
DEcomb$FC.HE.HU <- (DEcomb$Lin.FC.HE/DEcomb$Lin.FC.HU)
DEcomb$ratioZscore <- ((DEcomb$FC.HE.HU - mean(DEcomb$FC.HE.HU))/sd(DEcomb$FC.HE.HU))

#Set thresholds for significance of differential expression by multiple comparisons (BY)
DE_sigBY <- subset(DEcomb, BY.pvalue.HU < 0.05)
DE_sigBY <- subset(DEcomb, BY.pvalue.HE < 0.05)

#Set limits for z scores to determine genes with different expression between groups
DE_altResUnpair <- subset(DE_sigBY, ratioZscore > 1 | ratioZscore < -1)
write.csv(DE_altRes, "Desktop/DE_altRes.csv", row.names = FALSE)