rm(list = ls())
setwd("~/Desktop/Dissertation/Publications/PLL5 paper/Manuscript/Source data/Code/")
getwd()
library(tidyr)
library(tidyverse)
library(ggplot2)
my_data <- readRDS("../Code/Data/PPase linear regression plot/PPase activity_biorep6.rds")
PPase_data <- pivot_longer(my_data,cols = colnames(my_data)[3]:colnames(my_data)[8],values_to = "Absorbance")
PPase_data$Sample <- factor(PPase_data$Sample,levels=c("WT","3A","3D","freeGFP"))
colnames(PPase_data)[1] <- "SlPLL2" #PPase

ggplot(PPase_data%>%filter(time!="0h"), aes(x = time, y = Absorbance,color = SlPLL2,group = SlPLL2))+#facet_wrap(vars(Sample),nrow = 1)+
  geom_point()+stat_smooth(method = "lm",formula = y ~ x,geom = "smooth",alpha = 0.1)+scale_color_manual(values=c("#7570b3", "#e7298a", "#1b9e77","#e6ab02"))+
  ylab(label = "Absorbance of inorganic Pi")+scale_y_continuous(limits = c(-0.05,0.5))+#xlab("Time")+labs(title ="PPase activity")+
  theme(plot.title = element_text(size=14,color = "black",hjust = 0,face = "bold"))+
  theme(title = element_text(size=13,color = "black",face = "bold"))+
  theme(axis.text = element_text(size=12,color = "black",face = "bold"))+
  theme(panel.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),panel.grid.major = element_blank(),
        legend.background = element_rect(fill = "transparent",color = NA),
        legend.key = element_blank(),legend.box.background = element_rect(fill = "transparent",color = NA),legend.text=element_text(size=10))
#ggsave("../Code/Figure/PPase activity assay_biorep6.png", width=4.5,height=3)

dev.off()

