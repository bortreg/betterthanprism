#Load packages
library(ggplot2)
library(dplyr)

#Open Data from Exported from FlowJo and format for analysis
tfh <- read.csv("~/Desktop/20190917_TfhPanel_comb_v02.csv")
tfh$macaque.ID <- as.factor(tfh$macaque.ID)

#This is a way to "find and replace" with dplyr! Better than xcel :-)
tfh <- tfh %>%
  mutate(study.wk = replace(study.wk, study.wk == 16, 16.43))

#Index tfh to identify populations of interest
cd8Cxcr5 <- tfh[tfh$population == "Cells/Singlets/Live/Tcells/CD4-,CD8+/CXCR5+", ]
cd4pd1Dneg <- tfh[tfh$population == "Cells/Singlets/Live/Tcells/CD4+,CD8-/PD-1low/CXCR3-,CXCR5-", ]
cd4pd1X5 <- tfh[tfh$population == "Cells/Singlets/Live/Tcells/CD4+,CD8-/PD-1low/CXCR3-,CXCR5+", ]
cd4pd1X3 <- tfh[tfh$population == "Cells/Singlets/Live/Tcells/CD4+,CD8-/PD-1low/CXCR3+,CXCR5-", ]
cd4pd1Dpos <- tfh[tfh$population == "Cells/Singlets/Live/Tcells/CD4+,CD8-/PD-1low/CXCR3+,CXCR5+", ]

#Index df with specific timepoint
wk2 <- cd4pd1X5[cd4pd1X5$study.wk == -2, ]

#plot data for each population at a given timepoint
p <- ggplot(cd4pd1Dpos, aes(group = study.wk, x = study.wk, y = statistic)) +
  geom_boxplot(size = 0.3, width = 1.8, fill = "white") + theme_light()
wkP <- p + geom_jitter(width = 0.08, aes(shape = vaccine.route, color = macaque.ID))
wkPl <- wkP + geom_vline(xintercept = 0, linetype="solid", color = "black", size=0.3)+ 
  geom_vline(xintercept=8, linetype="dashed", color = "black", size = 0.3)+
  geom_vline(xintercept=16, linetype="dashed", color = "black", size = 0.3)

#Save PDF of plot
#Add plot labels
wkPlab <- wkPl + labs(title = "Circulating Tfh",subtitle = "PBMC", color = "macaque",
                   x = "Study Week", y = "%CXCR3+ CXCR5+ of PD-1low CD4 T cells")
wkPlab  

#Asign variables for column statistics for timepoints
col1 <- cd4pd1Dpos[cd4pd1Dpos$study.wk == -2, ]
col2 <- cd4pd1Dpos[cd4pd1Dpos$study.wk == 2, ]

#Check normality of data (data is normal if p > 0.05)
shapiro.test(col1$statistic)
qqnorm(col1$statistic,main = "Column 1", pch=19)
qqline(col1$statistic)

shapiro.test(col2$statistic)
qqnorm(col2$statistic,main = "Column 2", pch=19)
qqline(col2$statistic)




