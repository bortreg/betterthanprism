##Setup column graphs and run basic stats. This script is written for analyzing longitudinal data
##from NHP studies in which data is collected from microscopy, flow cytometry, ELISAs etc.

library(ggplot2)
library(dplyr)

#Import Data arranged with nhp id as rows and experiment and timepoint in columns

ovs3db <- read.csv("Desktop/OVS3db_v01.csv")

#Fix vector data types
ovs3db$mac.id <- as.character(ovs3db$mac.id)
ovs3db$timepoint <- as.character(ovs3db$timepoint)
ovs3db[,4:ncol(ovs3db)] <- 

#Make column for groups 
ovs3db$group <- c(rep_len("IM",6),(rep_len("IEP",6)))
ovs3db <- ovs3db %>% relocate(group, .after = mac.id) 

Add additional stats columns
ovs3db <- ovs3db %>% mutate()

#Set population of interest
poi <- ovs3db$per.tfh.cd4.flo.drln

#plot data in ggplot2
ggplot(ovs3db,aes(timepoint, poi))+
  geom_jitter(width = 0.2, size = 0.75, aes(color = group))+
  geom_errorbar(mapping = aes(x=timepoint, ymin=(poi + SD)), ymax=poi - SD)+
  geom_line(size = 0.4)+
  facet_wrap(~group)+
  theme_bw()
