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

ggplot(HU_DE, aes(x=Log2.fold.change, y=-log2(P.value), color = FDR)) +
  geom_point() +
  scale_color_manual(values = c("#45A242","#5CA25A","#89BE88","#C6E0C5","#CECFCE"))



#Merge dataframes
DEcomb <- merge(HU_DE, HEU_DE)

#Calculate z scores from ratios of fold change from experimental group over control group
DEcomb$FC.HE.HU <- (DEcomb$Lin.FC.HE/DEcomb$Lin.FC.HU)
DEcomb$ratioZscore <- ((DEcomb$FC.HE.HU - mean(DEcomb$FC.HE.HU))/sd(DEcomb$FC.HE.HU))

#Set thresholds for significance of differential expression by multiple comparisons (BY)
DE_sigBY <- subset(DEcomb, BY.pvalue.HU < 0.05)
DE_sigBY <- subset(DEcomb, BY.pvalue.HE < 0.05)

#Set limits for z scores to determine genes with different expression between groups
DE_altResUnpair <- subset(DE_sigBY, ratioZscore > 1 | ratioZscore < -1)
write.csv(DE_altRes, "Desktop/DE_altRes.csv", row.names = FALSE)