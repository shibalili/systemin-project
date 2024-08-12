rm(list = ls())
setwd("~/Desktop/Dissertation/Publications/PLL5 paper/Manuscript/Source data/Code/")
getwd()
library(tidyverse)
library(ggplot2)
Order_a <- readRDS("../Code/Data/k5_imp10_Order_a.rds")
ppep.y1 <- readRDS("../Code/Data/ORA/imp10_allBin1_k5-HA.rds")
table(ppep.y1%>%filter(pvalue<0.05)%>%select(sig_cluster))
ppep.y1.1 <- readRDS("../Code/Data/ORA/imp10_allBin1.1_k5-HA.rds")%>%filter(pvalue<0.05,sig_cluster=="k1"|sig_cluster=="k4")
df1.1 <- merge(Order_a[Order_a$ppep.site%in%ppep.y1.1$ppep.site,],ppep.y1.1,by ="ppep.site")%>%arrange(sig_cluster)%>%filter(pvalue<0.05,genotype=="WT")
df1.1$bin.name[str_detect(df1.1$bin.name, pattern = "20")] <- "bin 20 Cytoskeleton organisation"  
count1.1 <- plyr::count(df1.1, vars =c("bin.name","genotype","cluster","pvalue"))%>%
  mutate(bin.code=word(bin.name,1,2,sep = fixed(" ")))%>%arrange(cluster)

anno <- readRDS("../Code/Data/ORA/BackID2155_Annotation_detail_bin1.1.rds")
df <- merge(Order_a[Order_a$ppep.site%in%ppep.y1$ppep.site,],ppep.y1,by ="ppep.site")%>%arrange(sig_cluster)%>%filter(pvalue<0.05,genotype=="WT")
colnames(df)[12] <- "bin.name_short"
df <- df%>%mutate(bin.name=df$bin.name_short,.before = 1)

dd <- df%>%filter_at(vars(bin.name_short),any_vars(str_detect(., pattern = "18")))
bincode <- data.frame(Protein.ID=unique(dd$Protein.ID),bincode=anno[unique(dd$Protein.ID),3])
for (i in 1:nrow(bincode)) {bincode$bincode[i] <- unlist(strsplit(bincode[i,2],split=';', fixed=TRUE))[str_detect(strsplit(bincode[i,2],split=';', fixed=TRUE)%>%unlist(), pattern = "18")]%>%str_trim() }
bincode$bincode_short <- word(bincode$bincode,1,2,sep = fixed("."))
name <- readRDS("../Code/Data/ORA/Annotation source/ITAG3.2 mapman_Ensambl44.rds")%>%select(2:3)%>% filter(BINCODE %in% bincode$bincode_short)
colnames(name)[1] <- "bincode_short"
bincode <- merge(bincode,name,by = "bincode_short")
bincode$bin.name <- paste("bin",bincode$bincode_short,bincode$NAME)
for (j in 1:nrow(bincode)) {df[df$Protein.ID==bincode$Protein.ID[j],"bin.name"] <- bincode[j,"bin.name"]}

dd <- df%>%filter_at(vars(bin.name_short),any_vars(str_detect(., pattern = "24")))
bincode <- data.frame(Protein.ID=unique(dd$Protein.ID),bincode=anno[unique(dd$Protein.ID),3])
for (i in 1:nrow(bincode)) {bincode$bincode[i] <- unlist(strsplit(bincode[i,2],split=';', fixed=TRUE))[str_detect(strsplit(bincode[i,2],split=';', fixed=TRUE)%>%unlist(), pattern = "24")]%>%str_trim() }
bincode$bincode_short <- word(bincode$bincode,1,2,sep = fixed("."))
name <- readRDS("../Code/Data/ORA/Annotation source/ITAG3.2 mapman_Ensambl44.rds")%>%select(2:3)%>% filter(BINCODE %in% bincode$bincode_short)
colnames(name)[1] <- "bincode_short"
bincode <- merge(bincode,name,by = "bincode_short")
bincode$bin.name <- paste("bin",bincode$bincode_short,bincode$NAME)
for (j in 1:nrow(bincode)) {df[df$Protein.ID==bincode$Protein.ID[j],"bin.name"] <- bincode[j,"bin.name"]}
count <- rbind(plyr::count(df, vars =c("bin.name","genotype","cluster","pvalue"))%>%
                 mutate(bin.code=word(bin.name,1,2,sep = fixed(" ")))%>%arrange(cluster),count1.1)%>%arrange(cluster)
count$bin.name[grep("bin 18.10",count$bin.name)]<-"bin 18.10 Protein modification.peptide maturation"
count$bin.name[grep("bin 19.3",count$bin.name)]<- "bin 19.3 Protein autophagy"
count$bin.name[grep("bin 11.1",count$bin.name)]<- "bin 11.1 Phytohormone.ABA"
count$bin.name[grep("bin 24.1",count$bin.name)]<-"bin 24.1 Solute primary active transport"
count$bin.name[grep("bin 24.2",count$bin.name)]<-"bin 24.2 Solute carrier-mediated transport"
count$bin.name[grep("bin 27.3",count$bin.name)]<-"bin 27.3 SnRK1-kinase regulatory system"

count$bin.name <- factor(count$bin.name,levels=unique(count$bin.name)) # order bin name
ggplot(count,aes(cluster,bin.name))+geom_point(aes(color=pvalue,size=freq),show.legend = T)+scale_y_discrete(limits=rev)+
  #scale_x_discrete(labels=c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5"))+#xlim("1","2","3","4","5")+facet_grid(~as.factor(genotype))+
  scale_color_gradientn(name ="P value",colours=c("#371ea3","#46bac2", "#b3eebe"),guide=guide_colorbar(reverse=T))+
  scale_size_area(max_size = 8)+ #guides(size="none")+
  labs(y='',x='Phosphorylation clusters',size='P_sites number')+
  theme(axis.title=element_text(size=10,colour = 'black'), axis.text=element_text(size=10,colour = 'black'), 
        panel.background = element_rect(color='black',linewidth = 0.5,fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major.y = element_line(colour = "gray",linetype = "dotted"),
        panel.grid.major.x = element_blank(),
        legend.key = element_blank(),legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent",color = NA)) # get rid of legend panel bg
#ggsave("../Code/Figure/Github/ORA_imp10_K5_mix_pvalue.png",width=6.5,height =4)

dev.off()