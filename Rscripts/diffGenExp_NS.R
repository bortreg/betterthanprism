#Load Packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(calibrate)
library(RColorBrewer)

#Import Tables of Differentially Expressed Genes from nSolver
#*Make sure to rename columns you wish to remain unique to a specific group
HU_DE <- read.csv("~/Documents/R/Nanostring_DataFrames/DE_HU.csv")
HEU_DE <- read.csv("~/Documents/R/Nanostring_DataFrames/DE_HEU.csv")

#Plot heatmaps
df <- HEU_DE
tiff("~/Desktop/test.tiff", units = "in", width=5, height=5, res = 300)
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
DEtop100 <- DEcomb[1:100, c("HU.Log2.fold.change","HEU.Log2.fold.change")] 
row.names(DEtop100) <- DEcomb[1:100,"Gene"]
DEtop100 <- t(DEtop100)  

#Heatmap of DEtop100
zscore <- function(x) {          #zscore function
  z <- (mat_data - mean(x)) / sd(x)
  return(z)
}



#Calculate z scores from ratios of fold change from experimental group over control group
DEcomb$FC.HE.HU <- (DEcomb$Lin.FC.HE/DEcomb$Lin.FC.HU)
DEcomb$ratioZscore <- ((DEcomb$FC.HE.HU - mean(DEcomb$FC.HE.HU))/sd(DEcomb$FC.HE.HU))

#Set thresholds for significance of differential expression by multiple comparisons (BY)
DE_sigBY <- subset(DEcomb, BY.pvalue.HU < 0.05)
DE_sigBY <- subset(DEcomb, BY.pvalue.HE < 0.05)

#Set limits for z scores to determine genes with different expression between groups
DE_altResUnpair <- subset(DE_sigBY, ratioZscore > 1 | ratioZscore < -1)
write.csv(DE_altRes, "Desktop/DE_altRes.csv", row.names = FALSE)