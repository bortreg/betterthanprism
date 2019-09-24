#Load packages
library(ggplot2)
library(dplyr)


#Open Data from Exported from FlowJo
rawInnFlow <- read.csv("~/Desktop/ovs2FlowData/OVS2_innateFlow_Combined.csv")

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
lgLab <- lg + labs(title = "B cells",subtitle = "PBMC", color = "macaque",
          x = "Days post-vaccine", y = "Frequency of mDC") +
pdf("Desktop/graphBucket/OVS2flowPlots/DC_CD83_PBMC.pdf") +
lgLab +
dev.off()






#Graph a test population
DC83_InnFlow <- rawInnFlow[rawInnFlow$cell.population == "Cells/Singlets/HLADR+CD3-/Bcells/CD11c+", ]
DC83_PBMC <- DC83_InnFlow[DC83_InnFlow$tissue == "PBMC", ]
ggplot(DC83_PBMC, aes(x = study.day, y = statistic, color = macaque.ID, group = macaque.ID)) + 
  geom_point() + geom_line()

#Graph Multiple Tissues at once
lg <- ggplot(DC83_InnFlow, aes(x = study.day, y = statistic, 
                         color =  macaque.ID, 
                         shape =  vacc.route,
                         linetype = vacc.route,
                         group = macaque.ID)) + 
  scale_shape_manual(name = "Vaccine Route", values = c(19,1)) +
  geom_point() + geom_line() +
  facet_grid(~tissue, scales = "free_x", space = "free_x") +
  geom_vline(xintercept = 0, linetype="solid", color = "black", size=0.3)+ 
  geom_vline(xintercept=56, linetype="dashed", color = "black", size = 0.3)+
  geom_vline(xintercept=112, linetype="dashed", color = "black", size = 0.3)

#Add plot labels
lg + labs(title = "CD83+ mDC", color = "macaque")




