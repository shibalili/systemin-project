rm(list = ls())
setwd("~/Desktop/Dissertation/Publications/PLL5 paper/Manuscript/Source data/Code/")
getwd()
library(tidyverse)
library(ggplot2)
library(limma)

## Detect significant phosphorylation changes of phosphosites from time point 1-45 min against time 0

int_norm <- readRDS("../Code/Data/int_norm.rds")%>%select(c(2:37)) #WT
setup <- readRDS("../Code/Data/setup_WTsyr1.rds")%>%filter(genotype=="WT")
#int_norm <- readRDS("../Code/Data/int_norm.rds")%>%select(c(38:73)) #syr1
#setup <- readRDS("../Code/Data/setup_WTsyr1.rds")%>%filter(genotype!="WT")
levels(setup$time.point)# time 0 as reference: Compare to time 0
design <- model.matrix( ~ time.point, data = setup)
colSums(design)
fit <- lmFit(int_norm, design)#fit a lm
fit <- eBayes(fit)#t-test and so on
res1 <- topTable(fit, coef = "time.point1", number = Inf,adjust.method = "BH")
res2 <- topTable(fit, coef = "time.point2", number = Inf,adjust.method = "BH")
res5 <- topTable(fit, coef = "time.point5", number = Inf,adjust.method = "BH")
res15 <- topTable(fit, coef ="time.point15", number = Inf,adjust.method = "BH")
res45 <- topTable(fit, coef ="time.point45", number = Inf,adjust.method = "BH")
head(res1)
test_df <- rbind(data.frame(res1,P_site_ID=rownames(res1),Sig_time_pointX="1"),
                 data.frame(res2,P_site_ID=rownames(res2),Sig_time_pointX="2"),
                 data.frame(res5,P_site_ID=rownames(res5),Sig_time_pointX="5"),
                 data.frame(res15,P_site_ID=rownames(res15),Sig_time_pointX="15"),
                 data.frame(res45,P_site_ID=rownames(res45),Sig_time_pointX="45"))
test_df$Sig_time_pointX <- factor(test_df$Sig_time_pointX, levels =c("1","2","5","15","45"))
test_df$genotype <- "WT"
head(test_df)
rownames(test_df) <- NULL
time1 <- na.omit(res1[res1[,"P.Value"]<0.05,])
time1$P_site_ID <- rownames(time1)
time1$Sig_time_pointX <- "1"
time1$regulation <- ifelse(time1$logFC>0,"Up","Down")# P and deP Trend determine by LogFC

time2 <- na.omit(res2[res2[,"P.Value"]<0.05,])
time2$P_site_ID <- rownames(time2)
time2$Sig_time_pointX <- "2"
time2$regulation <- ifelse(time2$logFC>0,"Up","Down")

time5 <- na.omit(res5[res5[,"P.Value"]<0.05,])
time5$P_site_ID <- rownames(time5)
time5$Sig_time_pointX <- "5"
time5$regulation <- ifelse(time5$logFC>0,"Up","Down")

time15 <- na.omit(res15[res15[,"P.Value"]<0.05,])
time15$P_site_ID <- rownames(time15)
time15$Sig_time_pointX <- "15"
time15$regulation <- ifelse(time15$logFC>0,"Up","Down")

time45 <- na.omit(res45[res45[,"P.Value"]<0.05,])
time45$P_site_ID <- rownames(time45)
time45$Sig_time_pointX <- "45"
time45$regulation <- ifelse(time45$logFC>0,"Up","Down")

head(time1)
df <-rbind(time1[,c(4:5,7:9)],time2[,c(4:5,7:9)],time5[,c(4:5,7:9)],time15[,c(4:5,7:9)],time45[,c(4:5,7:9)])  
df$Sig_time_pointX <- factor(df$Sig_time_pointX, levels =c("1","2","5","15","45"))
table(df$Sig_time_pointX)
df$genotype <- "WT" #"syr1"
head(df)
rownames(df) <- NULL
df <- df[,c(3:6,1:2)]
table(df$regulation)
#saveRDS("../Code/Data/WT_pvalue.rds",object = df)

df$genotype <- "syr1"
#saveRDS("../Code/Data/syr1_pvalue.rds",object = df)
