#Load packages
library(ggplot2)

#Open Data from Exported from FlowJo
rawInnFlow <- read.csv("/Volumes/sodora_d/sodora/(08) Personal Folders/Matthew W/OVS/OVS2_DCmonocyte_v02_KF.csv")

#Reformat Factors and Levels for analysis
rawInnFlow$macaque.ID <- as.factor(rawInnFlow$macaque.ID)
levels(rawInnFlow$tissue) <- c("LNax", "LNax", "LNsm", "PBMC")
rawInnFlow$collection.date <- as.Date(rawInnFlow$collection.date, tryFormats = c("%Y%m%d"),
                                      optional = FALSE)

#Gave up on trying to make study.days column based on study.wk. Sadly continue in excel.
write.csv(rawInnFlow, "Desktop/rawInnFlow_1.csv")
rawInnFlow <- read.csv("Desktop/rawInnFlow_1.csv")
rawInnFlow$macaque.ID <- as.factor(rawInnFlow$macaque.ID)

#Identify a PBMC population
Bcells <- rawInnFlow[rawInnFlow$cell.population == "B cells", ]

#Graph PBMC
testGraph <- ggplot(Bcells, aes(x = study.day, y = statistic, 
                               color =  macaque.ID, 
                               shape =  vacc.route,
                               linetype = vacc.route,
                               group = macaque.ID)) + 
  scale_shape_manual(name = "Vaccine Route", values = c(19,1)) +
  geom_point() + geom_line() +
  geom_vline(xintercept = 0, linetype="solid", color = "black", size=0.3)+ 
  geom_vline(xintercept=56, linetype="dashed", color = "black", size = 0.3)+
  geom_vline(xintercept=112, linetype="dashed", color = "black", size = 0.3)

#Add plot labels
lgLab <- lg + labs(title = "B cells",subtitle = "PBMC", color = "macaque",
          x = "Days post-vaccine", y = "Frequency of mDC") +
pdf("Desktop/graphBucket/OVS2flowPlots/DC_CD83_PBMC.pdf") +
lgLab +
dev.off()






#Graph a test population
DC83_InnFlow <- rawInnFlow[rawInnFlow$cell.population == "DC/CD11c+ mDC/CD83+", ]
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




