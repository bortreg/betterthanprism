#Load packages
library(ggplot2)
library(dplyr)

#Open Data from Exported from FlowJo and format for analysis
tfh <- read.csv("~/Documents/RawData/OVS/20191015_TfhPanel_comb_v03.csv")
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
cd4pd1X5_I <- tfh[tfh$population == "Cells/Singlets/Live/Tcells/CD4+,CD8-/PD-1low/CXCR3-,CXCR5+/ICOS+", ]
cd4pd1X5_M <- tfh[tfh$population == "Cells/Singlets/Live/Tcells/CD4+,CD8-/PD-1low/CXCR3-,CXCR5+/Maf+", ]
cd4pd1low <- tfh[tfh$population == "Cells/Singlets/Live/Tcells/CD4+,CD8-/PD-1low", ]
cd4all <- tfh[tfh$population == "Cells/Singlets/Live/Tcells/CD4+,CD8-", ]

 #Sort data by macaque ID for paired comparisons
cd8Cxcr5 <- cd8Cxcr5[order(cd8Cxcr5$macaque.ID), ]
cd4pd1Dneg <- cd4pd1Dneg[order(cd4pd1Dneg$macaque.ID), ]
cd4pd1X5 <- cd4pd1X5[order(cd4pd1X5$macaque.ID), ]
cd4pd1X3 <- cd4pd1X3[order(cd4pd1X3$macaque.ID), ]
cd4pd1Dpos <- cd4pd1Dpos[order(cd4pd1Dpos$macaque.ID), ]
cd4pd1X5_M <- cd4pd1X5_M[order(cd4pd1X5_M$macaque.ID), ]
cd4pd1X5_I <- cd4pd1X5_I[order(cd4pd1X5_I$macaque.ID), ]
cd4pd1low <- cd4pd1low[order(cd4pd1low$macaque.ID), ]
cd4all <- cd4all[order(cd4all$macaque.ID), ]

#determine averages for for populations of interest by group
PoI <- cd4pd1X5_I
PoI_IM <- PoI[PoI$vaccine.route == "IM", ]
PoI_IEP <- PoI[PoI$vaccine.route == "Oral", ]

#calculate percent circulating tfh of total CD4 T cells
cd4all$statTFH <- ((cd4pd1Dpos$cell.num + cd4pd1X5$cell.num)/cd4all$cell.num)*100

#plot data for each population at a given timepoint
p <- ggplot(PoI, aes(x = study.wk, y = statistic, color = macaque.ID))
p + geom_line(size = 0.4) + 
  theme_bw() +
  stat_smooth(aes(group = vaccine.route)) +
  stat_summary(aes(group = vaccine.route), geom = "point", fun = mean, shape = 16, size = 1) +
  facet_grid(. ~vaccine.route)

#Add plot labels
wkPlab <- wkPl + labs(title = "Circulating Tfh",subtitle = "PBMC", color = "macaque",
                   x = "Study Week", y = "%ICOS+ of CXCR3- cTfh cells")
wkPlab  


#Asign variables for column statistics for timepoints
col1 <- cd4pd1X5_I[cd4pd1X5_I$study.wk == -2, ]
col1oral <- col1[col1$vaccine.route == "Oral", ]
col1im <- col1[col1$vaccine.route == "IM", ]
col2 <- cd4pd1X5_I[cd4pd1X5_I$study.wk == 16.43, ]
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


