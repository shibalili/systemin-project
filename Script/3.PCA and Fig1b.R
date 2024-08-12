rm(list = ls())
setwd("~/Desktop/Dissertation/Publications/PLL5 paper/Manuscript/Source data/Code/")
getwd()
library(tidyverse)
library(ggplot2)
###### Ave prep ######
#WT sys
int_norm <- readRDS("../Code/Data/int_norm.rds")%>%select(c(2:37)) 
WT_S0 <- apply(int_norm[,1:6], 1, function(x){return(mean(x, na.rm = TRUE))})
WT_S1 <- apply(int_norm[,7:12], 1, function(x){return(mean(x, na.rm = TRUE))})
WT_S2 <- apply(int_norm[,13:18], 1, function(x){return(mean(x, na.rm = TRUE))})
WT_S5 <- apply(int_norm[,19:24], 1, function(x){return(mean(x, na.rm = TRUE))})
WT_S15 <-apply(int_norm[,25:30], 1, function(x){return(mean(x, na.rm = TRUE))})
WT_S45 <-apply(int_norm[,31:36], 1, function(x){return(mean(x, na.rm = TRUE))})
Smean_WT <- cbind(WT_S0,WT_S1,WT_S2,WT_S5,WT_S15,WT_S45)%>%as.data.frame()

#syr1 sys
int_norm <- readRDS("../Code/Data/int_norm.rds")%>%select(c(38:73))
syr1_S0 <- apply(int_norm[,1:6], 1, function(x)mean(x, na.rm = TRUE))
syr1_S1 <- apply(int_norm[,7:12], 1, function(x)mean(x, na.rm = TRUE))
syr1_S2 <- apply(int_norm[,13:18], 1,function(x) mean(x, na.rm = TRUE))
syr1_S5 <- apply(int_norm[,19:24], 1, function(x)mean(x, na.rm = TRUE))
syr1_S15 <- apply(int_norm[,25:30], 1, function(x)mean(x, na.rm = TRUE))
syr1_S45 <- apply(int_norm[,31:36], 1, function(x)mean(x, na.rm = TRUE))
Smean_syr1 <- cbind(syr1_S0,syr1_S1,syr1_S2,syr1_S5,syr1_S15,syr1_S45)%>%as.data.frame()
#identical(rownames(Smean_WT),rownames(Smean_syr1))
Smean <- cbind(Smean_WT,Smean_syr1)
#saveRDS("../Code/Data/int_Smean.rds",object = Smean)

## remove P_sites with same significant phosphorylation in both WT and syr1
WT <- readRDS("../Code/Data/WT_pvalue.rds")
WT$combine <- paste(WT$P_site_ID,WT$regulation,WT$Sig_time_pointX)
table(WT$Sig_time_pointX)
length(unique(WT$P_site_ID))#911
syr1 <- readRDS("../Code/Data/syr1_pvalue.rds")
syr1$combine <- paste(syr1$P_site_ID,syr1$regulation,syr1$Sig_time_pointX)
WT%>%filter(combine%in%intersect(WT$combine, syr1$combine))%>%dim() #remove shared DEGs from WT sides: -61
p <- NULL
pp <- NULL
ppp <- NULL
a <- NULL
for (j in 1:length(levels(WT$Sig_time_pointX))) {
  p <- WT%>%filter(Sig_time_pointX==levels(WT$Sig_time_pointX)[j])
  pp <-syr1$combine
  a <- p[p$combine %in% pp,]
  ppp <- rbind(ppp,p[!p$combine%in%a$combine,])
}
table(ppp$Sig_time_pointX)
length(unique(ppp$P_site_ID))#883
#saveRDS("../Code/Data/DEG_WTremsyr1.rds",object = ppp)#1299

##### PCA #####
setup <- readRDS("../Code/Data/setup_Smean.rds")
Smean <- readRDS("../Code/Data/int_Smean.rds")
ppp <- readRDS("../Code/Data/DEG_WTremsyr1.rds")
Smean <- Smean[intersect(rownames(Smean),unique(ppp$P_site_ID)),]
pca <- prcomp(x = t(na.omit(Smean)),center =TRUE,scale=TRUE)
summary(pca)
df_pca<- as.data.frame(pca$x)
df_pca$time.point <- setup[rownames(setup),2]
df_pca$gentoype <- setup[rownames(df_pca),3]

ggplot(df_pca,aes(x=PC1,y=PC2,color =gentoype,label=row.names(df_pca),shape=time.point))+
  geom_point(size=3,show.legend = T, stroke=1)+scale_colour_manual(name="Genotype",values= c("#dc1400", "#00AFBB"))+
  ggtitle("PC1 and PC2",subtitle = "systemin:WT vs syr1")+
  scale_shape_manual(name="Time (min)",values=c(15,16,17,0,1,2))+
  scale_y_continuous(limits=c(-12, 14))+scale_x_continuous(limits=c(-20, 18))+
  xlab(label = paste0("PC1:",round(summary(pca)[[6]][2,1]*100,2),"% variance"))+
  ylab(label = paste0("PC2:",round(summary(pca)[[6]][2,2]*100,2),"% variance"))+
  guides(color = guide_legend(order=1),shape = guide_legend(order=2))+
  theme(title = element_text(size=16,color = "black"),
        axis.text = element_text(size=14,color = "black"),axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        legend.position="right",
        legend.background = element_rect(fill = "transparent",color = NA),legend.key = element_blank(),
        legend.box.background = element_rect(fill = "transparent",color = NA),legend.text=element_text(size=14))
ggsave(filename = "../Code/Figure/Github/PCA_DEG_PC12.png",width =5,height =4)
dev.off()
