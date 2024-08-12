rm(list = ls())
setwd("~/Desktop/Dissertation/Publications/PLL5 paper/Manuscript/Source data/Code/")
getwd()
library(tidyverse)
## clean data
dat <- read.table("../Code/Data/MQ_txt/Phospho (STY)Sites.txt",header=TRUE,sep = "\t",quote="",fill = FALSE)#6179
dat <- dat%>%mutate(ppep.site=paste(substr(dat[,4],1,18),dat[,'Positions.within.proteins']),.before=1)
dat <- dat[-grep("+",dat[,"Reverse"], fixed=TRUE),]#6123
dat <- dat[-grep("+",dat[,"Potential.contaminant"], fixed=TRUE),]#6112
dat <- dat[dat[,"Localization.prob"]>0.75,]
rownames(dat) <- dat$ppep.site
#saveRDS("../Code/Data/dat_cleaned.rds",object = dat)

## Normalize and log2 intensity
dat <- readRDS("../Code/Data/dat_cleaned.rds")
grep("Intensity.syr1_S01",colnames(dat))
colnames(dat)[grep("Intensity.syr1_S01",colnames(dat))] #389-Intensity.syr1_S01
grep("Intensity.WT_S56",colnames(dat))
colnames(dat)[grep("Intensity.WT_S56",colnames(dat))] #460-Intensity.WT_S56
int <- dat[,c(389:460)]
colnames(int)
setup <- readRDS("../Code/Data/setup_WTsyr1.rds")
int <- int[,rownames(setup)]
colnames(int)
int[int==0] <- NA
normalize.intensities <- function(x){
  x0 <- matrix(NA, ncol = ncol(x), nrow = nrow(x))
  overall.mean <- mean(apply(x, 2, function(x){return(sum(x, na.rm = TRUE))}))#average the sum intensity of each column/sample
  print(overall.mean)
  for(i in 1:ncol(x0)){
    x0[,i] <- overall.mean*x[,i]/sum(x[,i], na.rm = TRUE)
  }
  colnames(x0) <- colnames(x)
  rownames(x0) <- rownames(x)
  return(x0)
}
intensities<- normalize.intensities(int)%>%as.data.frame()%>%log2()%>%mutate(ppep.site=row.names(.),.before = 1)
dim(intensities)#4804 73
#saveRDS(file = "../Code/Data/int_norm.rds", object = intensities)


