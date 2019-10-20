#Load packages
library(ggplot2)
library(dplyr)

#Open Data from Exported from FlowJo and format for analysis
tfh <- read.csv("~/Desktop/20191015_TfhPanel_comb_v03.csv")
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

#plot data for each population at a given timepoint
p <- ggplot(cd4pd1X5, aes(group = study.wk, x = study.wk, y = statistic)) +
  geom_boxplot(size = 0.6, width = 0.6, fill = "lightgray") + theme_light()
wkP <- p + geom_jitter(width = 0.08, aes(shape = vaccine.route, color = macaque.ID))
wkPl <- wkP + geom_vline(xintercept = 0, linetype="solid", color = "black", size=0.3)+ 
  geom_vline(xintercept=8, linetype="dashed", color = "black", size = 0.3)+
  geom_vline(xintercept=16, linetype="dashed", color = "black", size = 0.3)

wkPl

#plot 2 time points for a given population
wkA <- cd4pd1Dpos[cd4pd1Dpos$study.wk == 8, ]
wkB <- cd4pd1Dpos[cd4pd1Dpos$study.wk == 8.43, ] 
p <- ggplot(wkA, aes(group = study.wk, x = study.wk, y = statistic)) +
  geom_boxplot(size = 0.3, width = 1.8, fill = "white") + theme_light()
wkP <- p + geom_jitter(width = 0.08, aes(shape = vaccine.route, color = macaque.ID))
wkPl <- wkP + geom_vline(xintercept = 0, linetype="solid", color = "black", size=0.3)

#Add plot labels
wkPlab <- wkPl + labs(title = "Circulating Tfh",subtitle = "PBMC", color = "macaque",
                   x = "Study Week", y = "%CXCR3- CXCR5+ of PD-1low CD4 T cells")
wkPlab  

#Sort data by macaque ID for paired comparisons
cd8Cxcr5 <- cd8Cxcr5[order(cd8Cxcr5$macaque.ID), ]
cd4pd1Dneg <- cd4pd1Dneg[order(cd4pd1Dneg$macaque.ID), ]
cd4pd1X5 <- cd4pd1X5[order(cd4pd1X5$macaque.ID), ]
cd4pd1X3 <- cd4pd1X3[order(cd4pd1X3$macaque.ID), ]
cd4pd1Dpos <- cd4pd1Dpos[order(cd4pd1Dpos$macaque.ID), ]


#Asign variables for column statistics for timepoints
col1 <- cd4pd1X5[cd4pd1X5$study.wk == 8, ]
col1oral <- col1[col1$vaccine.route == "Oral", ]
col1im <- col1[col1$vaccine.route == "IM", ]
col2 <- cd4pd1X5[cd4pd1X5$study.wk == 8.43, ]
col2oral <- col2[col2$vaccine.route == "Oral", ]
col2im <- col2[col2$vaccine.route == "IM", ]



#Check normality of data (data is normal if p > 0.05)
shapiro.test(col1$statistic)
qqnorm(col1$statistic,main = "Column 1", pch=19)
qqline(col1$statistic)
shapiro.test(col1oral$statistic)
qqnorm(col1oral$statistic,main = "Column 1", pch=19)
qqline(col1oral$statistic)
shapiro.test(col1im$statistic)
qqnorm(col1im$statistic,main = "Column 1", pch=19)
qqline(col1im$statistic)

shapiro.test(col2$statistic)
qqnorm(col2$statistic,main = "Column 2", pch=19)
qqline(col2$statistic)
shapiro.test(col2oral$statistic)
qqnorm(col2oral$statistic,main = "Column 1", pch=19)
qqline(col2oral$statistic)
shapiro.test(col2im$statistic)
qqnorm(col2im$statistic,main = "Column 1", pch=19)
qqline(col2im$statistic)

#ttest of data
t.test(col1$statistic, col2$statistic, paired = TRUE)
t.test(col1im$statistic, col2im$statistic, paired = TRUE)
t.test(col1oral$statistic, col2oral$statistic, paired = TRUE)


