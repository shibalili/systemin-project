rm(list = ls())
setwd("~/Desktop/Dissertation/Publications/PLL5 paper/Manuscript/Source data/Code/")
getwd()
library(tidyverse)
dat <- readRDS("../Code/Data/dat_cleaned.rds")
grep("Number.of.Phospho..STY.",colnames(dat)) #302
## sequence
seq <-dat%>%select(c(302:304,307))%>%mutate(Flanking = gsub("\\;.*$", "", Sequence.window)%>%str_sub(9,23)) 
colnames(seq)[1] <- "No.STY"
colnames(seq)[4] <- "Probability"
## normalized-log2 intensity
int_norm <- readRDS("../Code/Data/int_norm.rds")
## number of measurement
scount <- readRDS("../Code/Data/Scount.rds")
df <- bind_cols(seq,scount[-1],int_norm)%>%mutate(Protein.ID=substr(ppep.site,1,18),.after = ppep.site)%>%relocate(ppep.site,Protein.ID)
## annotation
mapping <- readRDS("../Code/Data/ORA/Annotation source/phytozome_ITAG3.2_standard.rds")
ss <- mapping[unique(df$Protein.ID),c(1:10)]
df <- merge(ss,df,by="Protein.ID")%>%relocate(ppep.site,Protein.ID)
rownames(df) <- df$ppep.site
df <- df[rownames(int_norm),] #order rownames as the same order of int_norm #4804 102
## DEG
DEG<- readRDS("../Code/Data/DEGplot_a.rds")%>%mutate(combine = paste(genotype,regulation,Sig_time_pointX))
table(DEG$genotype,DEG$genotype_share)
Order_DEG <- DEG[order(DEG$Sig_time_pointX),] 
Order_DEG <- DEG[order(DEG$genotype,decreasing = F),] # WT comes first
df <- df%>%mutate(DEG="no",.before = 1)
l <- list()
for (i in 1:nrow(df)) {
  l[[i]]<- Order_DEG[which(Order_DEG$P_site_ID==rownames(df)[i]),"combine"]
  names(l)[i] <- rownames(df)[i]}#identical(names(l),rownames(df))
df$DEG<- l%>%as.character()
df[df$DEG=="character(0)","DEG"] <- "no"
#the regex "^c\\(|\\)$" only affects the first/last chars in the string, remove "c(" and ")" from string
df$DEG <- df$DEG%>%gsub("\"", "", .)%>%gsub("^c\\(|\\)$", "",.) 
## Kmeans clustering
df <- df%>%mutate(syr1_cluster="filtered",.before = 1)
df <- df%>%mutate(WT_cluster="filtered",.before = 1)
distribution <- readRDS("../Code/Data/k5_imp10_dist.rds")#distribution <- distribution%>%filter(WT!=syr1)#remove WT_cluster==syr1_cluster
table(distribution$WT,distribution$syr1)
df[which(rownames(df) %in% rownames(distribution)),]$WT_cluster<- distribution[rownames(df[which(rownames(df) %in% rownames(distribution)),]),"WT"]
df[which(rownames(df) %in% rownames(distribution)),]$syr1_cluster<- distribution[rownames(df[which(rownames(df) %in% rownames(distribution)),]),"syr1"]
table(df$WT_cluster,df$syr1_cluster)
## annotation:X4.2_symbol
X4.2_symbol <- readRDS("../Code/Data/ORA/Annotation source/phytozome_ITAG3.2_X4.2.rds")
df <- df%>%mutate(X4.2_symbol=X4.2_symbol[df$Protein.ID,c("X4.2_symbol")])%>%relocate(X4.2_symbol,.before = arabi.defline)
## annotation:Localization
df$Ara_SUBA3 <- gsub("[\"]","",df$Ara_SUBA3)%>%noquote()
df$Localization <- "null"
df[which(str_detect(df$Ara_SUBA3, "plasma")),"Localization"] <- "PM"
df[which(str_detect(df$Ara_SUBA3, "cytosol")),"Localization"] <- "cytosol"
df[which(str_detect(df$Ara_SUBA3, "endoplasmic")),"Localization"] <- ifelse(df[which(str_detect(df$Ara_SUBA3, "endoplasmic")),"Localization"]=="null","ER",df[which(str_detect(df$Ara_SUBA3, "endoplasmic")),"Localization"])
df[which(str_detect(df$Ara_SUBA3, "nucleus")),"Localization"] <- ifelse(df[which(str_detect(df$Ara_SUBA3, "nucleus")),"Localization"]=="null","nucleus",df[which(str_detect(df$Ara_SUBA3, "nucleus")),"Localization"])
df[which(str_detect(df$Ara_SUBA3, "extracellular")),"Localization"] <- ifelse(df[which(str_detect(df$Ara_SUBA3, "extracellular")),"Localization"]=="null","extracellular",df[which(str_detect(df$Ara_SUBA3, "extracellular")),"Localization"])
df[which(str_detect(df$Ara_SUBA3, "plastid")),"Localization"] <- ifelse(df[which(str_detect(df$Ara_SUBA3, "plastid")),"Localization"]=="null","plastid",df[which(str_detect(df$Ara_SUBA3, "plastid")),"Localization"])
df[which(str_detect(df$Ara_SUBA3, "mitochondrion")),"Localization"] <- ifelse(df[which(str_detect(df$Ara_SUBA3, "mitochondrion")),"Localization"]=="null","mitochondrion",df[which(str_detect(df$Ara_SUBA3, "mitochondrion")),"Localization"])
df[which(str_detect(df$Ara_SUBA3, "golgi")),"Localization"] <- ifelse(df[which(str_detect(df$Ara_SUBA3, "golgi")),"Localization"]=="null","golgi",df[which(str_detect(df$Ara_SUBA3, "golgi")),"Localization"])
df[which(str_detect(df$Ara_SUBA3, "peroxisome")),"Localization"] <- ifelse(df[which(str_detect(df$Ara_SUBA3, "peroxisome")),"Localization"]=="null","peroxisome",df[which(str_detect(df$Ara_SUBA3, "peroxisome")),"Localization"])
df[which(str_detect(df$Ara_SUBA3, "vacuole")),"Localization"] <- ifelse(df[which(str_detect(df$Ara_SUBA3, "vacuole")),"Localization"]=="null","vacuole",df[which(str_detect(df$Ara_SUBA3, "vacuole")),"Localization"])
sort(table(df$Localization))
df<- df%>%relocate(Localization,.before = Ara_SUBA3)
## annotation:Protein.type
mapping <- readRDS("../Code/Data/ORA/Annotation source/phytozome_ITAG3.2_standard.rds")
a <- mapping[unique(df$Protein.ID),]
TF<- a%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., pattern = "transcription factor")))%>%
  unique()%>%pull(Protein.ID)
receptor <- a%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., pattern = "receptor")))%>%
  unique()%>%pull(Protein.ID)
kinase<- a%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., pattern = "kinase")))%>%
  unique()%>%pull(Protein.ID)
fake_kinase <- a%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., pattern = "kinase")))%>%
  filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., pattern = "phosphatase")))
##### more true kinase than PPase so PPase first ####
phostase<- a%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., pattern = "phosphatase")))%>%
  unique()%>%pull(Protein.ID)
df$Protein.type <- "unknown"
df[which(df$Protein.ID %in% TF),"Protein.type"] <- "TF"
df[which(df$Protein.ID %in% receptor),"Protein.type"] <- "receptor"
df[which(df$Protein.ID %in% phostase),"Protein.type"] <- "PPase"
df[which(df$Protein.ID %in% kinase),"Protein.type"] <- "kinase"
df[which(df$Protein.ID%in%c("Solyc03g118350.3.1","Solyc10g005640.3.1")),"Protein.type"] <- "PPase"
table(df$Protein.type)
df<- df%>%relocate(Protein.type,.after = Protein.ID)
colnames(df)
df <- df[,c(4:5,1:3,6:108)]
df <- df[,c(1:2,23:24,3:22,25:108)]
#saveRDS("../Code/Data/Summary_df.rds",object = df)
