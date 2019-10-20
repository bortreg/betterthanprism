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
poI <- rawInnFlow[rawInnFlow$population == "Cells/Singlets/HLADR+CD3-/CD20-/CD16-CD14-/mDC", ]
poI <- poI[poI$tissue == "LNsm", ]

#Graph Population of Interest
p <- ggplot(poI, aes(group = study.week, x = study.week, y = statistic)) +
  geom_boxplot(size = 0.6, width = 0.3, fill = "lightgray") + theme_light()
wkP <- p + geom_jitter(width = 0.08, aes(shape = vaccine.route, color = macaque.ID))
wkPl <- wkP + geom_vline(xintercept = 0, linetype="solid", color = "black", size=0.3) 
  

#Add plot labels
poIlab <- substring(as.character(poI[1,5]), 27, last = 100000L)
poIlab <- gsub(pattern = "/", replacement = "_", poIlab)
poIlabP <- paste("~/Desktop/graphBucket/OVS2_LNsm/", poIlab, ".pdf", sep = "") 
poIlabP <- gsub(pattern = "-", replacement = "neg", poIlabP)
poIlabP <- gsub(pattern = "+", replacement = "pos", poIlabP, fixed = TRUE)


wkPlab <- wkPl + labs(title = poIlab, subtitle = "LNsm", color = "macaque",
          x = "Weeks post-vaccine", y = "Frequency of parent")

pdf(file = poIlabP)
wkPlab
dev.off()





