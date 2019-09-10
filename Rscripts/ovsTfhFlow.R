#Load packages
library(ggplot2)

#Open Data from Exported from FlowJo and format for analysis
tfh <- read.csv("~/Desktop/20190905_TfhPanel_comb_v01.csv")
tfh$macaque.ID <- as.factor(tfh$macaque.ID)

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
  geom_boxplot(width = 1.8, fill = "white") + theme_light()
p + geom_jitter(width = 0.08, aes(color = macaque.ID, shape = vaccine.route))
