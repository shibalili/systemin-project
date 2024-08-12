rm(list = ls())
setwd("~/Desktop/Dissertation/Publications/PLL5 paper/Manuscript/Source data/Code/")
getwd()
library(tidyverse)
library(missForest)
library(ggplot2)

###### imputation ######## 
int_combine <- readRDS("../Code/Data/int_norm.rds")%>%select(2:73)
int <- t(int_combine[apply(int_combine, 1, function(x){sum(!is.na(x))>10}),]) #P_site has 10 counts or more
dim(int)#>10 2450 50%

set.seed(123)#in order to get reproducible results
int.imp <- missForest(xmis = int,verbose = TRUE)
str(int.imp)
#saveRDS("../Code/Data/int_norm_imp10.rds",object = int.imp) #missForest dataform 
int.imp <- readRDS("../Code/Data/int_norm_imp10.rds")
imp10 <- t(int.imp$ximp)%>%as.data.frame()
imp10$Protein.ID <- substr(rownames(imp10),1,18)
imp10$ppep.site <- rownames(imp10)
imp10 <- imp10[, c(74,73,1:72)]
#saveRDS("../Code/Data/imp10_df.rds",object = imp10)
###### imp10 Ave prep ######
#WT sys
int_norm <- readRDS("../Code/Data/imp10_df.rds")%>%select(c(3:38)) 
WT_S0 <- apply(int_norm[,1:6], 1, function(x){return(mean(x, na.rm = TRUE))})
WT_S1 <- apply(int_norm[,7:12], 1, function(x){return(mean(x, na.rm = TRUE))})
WT_S2 <- apply(int_norm[,13:18], 1, function(x){return(mean(x, na.rm = TRUE))})
WT_S5 <- apply(int_norm[,19:24], 1, function(x){return(mean(x, na.rm = TRUE))})
WT_S15 <-apply(int_norm[,25:30], 1, function(x){return(mean(x, na.rm = TRUE))})
WT_S45 <-apply(int_norm[,31:36], 1, function(x){return(mean(x, na.rm = TRUE))})
Smean_WT <- cbind(WT_S0,WT_S1,WT_S2,WT_S5,WT_S15,WT_S45)%>%as.data.frame()
#syr1 sys
int_norm <- readRDS("../Code/Data/imp10_df.rds")%>%select(c(39:74))
syr1_S0 <- apply(int_norm[,1:6], 1, function(x)mean(x, na.rm = TRUE))
syr1_S1 <- apply(int_norm[,7:12], 1, function(x)mean(x, na.rm = TRUE))
syr1_S2 <- apply(int_norm[,13:18], 1,function(x) mean(x, na.rm = TRUE))
syr1_S5 <- apply(int_norm[,19:24], 1, function(x)mean(x, na.rm = TRUE))
syr1_S15 <- apply(int_norm[,25:30], 1, function(x)mean(x, na.rm = TRUE))
syr1_S45 <- apply(int_norm[,31:36], 1, function(x)mean(x, na.rm = TRUE))
Smean_syr1 <- cbind(syr1_S0,syr1_S1,syr1_S2,syr1_S5,syr1_S15,syr1_S45)%>%as.data.frame()
Smean <- cbind(Smean_WT,Smean_syr1)
#saveRDS("../Code/Data/int_Smean_imp10.rds",object = Smean)
Smean_WT <- readRDS("../Code/Data/int_Smean_imp10.rds")%>%select(c(1:6))%>%mutate(genotype="WT")
rownames(Smean_WT) <- paste0(rownames(Smean_WT),"_",word(colnames(Smean_WT)[1],sep = fixed("_")))
colnames(Smean_WT)[1:6] <- word(colnames(Smean_WT)[1:6],2,sep = fixed("_"))
Smean_syr1 <- readRDS("../Code/Data/int_Smean_imp10.rds")%>%select(c(7:12))%>%mutate(genotype="syr1")
rownames(Smean_syr1) <- paste0(rownames(Smean_syr1),"_",word(colnames(Smean_syr1)[1],sep = fixed("_")))
colnames(Smean_syr1)[1:6] <- word(colnames(Smean_syr1)[1:6],2,sep = fixed("_"))
Smean_df <- na.omit(rbind(Smean_WT,Smean_syr1))
table(Smean_df$genotype)
#saveRDS("../Code/Data/kmc_imp10_df.rds",object = Smean_df)

