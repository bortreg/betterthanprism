#Load packages
library(ggplot2)
library(dplyr)


#Open Data from Exported from FlowJo
rawInnFlow <- read.csv("~/Desktop/ovsFlowData/OVS2/OVS2_innateFlow_Combined.csv")

#Reformat Factors and Levels for analysis
rawInnFlow$macaque.ID <- as.factor(rawInnFlow$macaque.ID)
levels(rawInnFlow$tissue) <- c("LNax", "LNing", "LNsm", "PBMC")
rawInnFlow$collection.date <- as.character(rawInnFlow$collection.date)
rawInnFlow$collection.date <- as.Date(rawInnFlow$collection.date, tryFormats = c("%Y%m%d"),
                                      optional = FALSE)

#Identify a population and tissue type
poI <- rawInnFlow[rawInnFlow$population == "Cells/Singlets/HLADR+CD3-/CD20-/CD16-CD14-/mDC/CD80+", ]
poI <- poI[poI$tissue == "PBMC", ]

#Graph Population of Interest
p <- ggplot(poI, aes(group = study.week, x = study.week, y = statistic)) +
  geom_boxplot(size = 0.3, width = 0.5, fill = "white") + theme_light()
wkP <- p + geom_jitter(width = 0.08, aes(shape = vaccine.route, color = macaque.ID))
wkPl <- wkP + geom_vline(xintercept = 0, linetype="solid", color = "black", size=0.3)+ 
  geom_vline(xintercept=8, linetype="dashed", color = "black", size = 0.3)+
  geom_vline(xintercept=16, linetype="dashed", color = "black", size = 0.3)


#Add plot labels
poIlab <- substring(as.character(poI[1,5]), 27, last = 100000L)

wkPlab <- wkPl + labs(title = poIlab, subtitle = "PBMC", color = "macaque",
          x = "Weeks post-vaccine", y = "Frequency of parent") +
pdf("Desktop/graphBucket/OVS2flowPlots/DC_CD83_PBMC.pdf") +
lgLab +
dev.off()






