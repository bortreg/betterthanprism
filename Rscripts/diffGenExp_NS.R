#Import Tables of Differentially Expressed Genes from nSolver
#*Make sure to rename columns you wish to remain unique to a specific group
HE_DE <- read.csv("Desktop/NS_analysisCBMC/HEunpair_DE.csv", sep = ",")
HU_DE <- read.csv("Desktop/NS_analysisCBMC/HUunpair_DE.csv", sep = ",")

#Merge dataframes
DEcomb <- merge(HE_DE, HU_DE)

#Calculate z scores from ratios of fold change from experimental group over control group
DEcomb$FC.HE.HU <- (DEcomb$Lin.FC.HE/DEcomb$Lin.FC.HU)
DEcomb$ratioZscore <- ((DEcomb$FC.HE.HU - mean(DEcomb$FC.HE.HU))/sd(DEcomb$FC.HE.HU))

#Set thresholds for significance of differential expression by multiple comparisons (BY)
DE_sigBY <- subset(DEcomb, BY.pvalue.HU < 0.05)
DE_sigBY <- subset(DEcomb, BY.pvalue.HE < 0.05)

#Set limits for z scores to determine genes with different expression between groups
DE_altResUnpair <- subset(DE_sigBY, ratioZscore > 1 | ratioZscore < -1)
write.csv(DE_altRes, "Desktop/DE_altRes.csv", row.names = FALSE)