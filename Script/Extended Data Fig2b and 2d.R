rm(list = ls())
setwd("~/Desktop/Dissertation/Publications/PLL5 paper/Manuscript/Source data/Code/")
getwd()
library(tidyverse)
library(ggplot2)
library(ggpubr) #get_legend
library(viridis)
#### PPase count #####
df <- readRDS("../Code/Data/Summary_df.rds")
table(df$Protein.type)
PPase <-df%>%filter(Protein.type=="PPase")%>%filter_at(vars(MapmanX4.2_Bincode),any_vars(str_detect(.,"18.")))
length(unique(PPase$Protein.ID)) #16
Count_ID_PPase <- plyr::count(PPase, vars =c("Protein.ID"))%>%arrange(freq)
PP2C<- df%>%filter(Protein.type=="PPase")%>%filter_at(vars(MapmanX4.2_Bincode),any_vars(str_detect(.,"18.")))%>%filter_at(vars("MapmanX4.2_Bin name",Ara_Mapman_DESCRIPTION),any_vars(str_detect(.,"PP2C")))%>%mutate(class="PPM/PP2C")
PPP <- df%>%filter(Protein.type=="PPase")%>%filter_at(vars(MapmanX4.2_Bincode),any_vars(str_detect(.,"18.")))%>%filter_at(vars("MapmanX4.2_Bin name",Ara_Mapman_DESCRIPTION),any_vars(str_detect(.,"PPP")))%>%mutate(class="PPP")
PTP <- df%>%filter(Protein.type=="PPase")%>%filter_at(vars(MapmanX4.2_Bincode),any_vars(str_detect(.,"18.")))%>%filter_at(vars("MapmanX4.2_Bin name",Ara_Mapman_DESCRIPTION),any_vars(str_detect(.,"PTP")))%>%mutate(class="PTP")
Asp <- df%>%filter(Protein.type=="PPase")%>%filter_at(vars(MapmanX4.2_Bincode),any_vars(str_detect(.,"18.")))%>%filter_at(vars("MapmanX4.2_Bin name",Ara_Mapman_DESCRIPTION),any_vars(str_detect(.,"aspartate")))%>%mutate(class="Asp-dep")
PPase_class <- rbind(PP2C[,c(1:2,109)],PPP[,c(1:2,109)],PTP[,c(1:2,109)],Asp[,c(1:2,109)])
count <- df%>%filter(ppep.site%in%PPase_class$ppep.site)

Count_PPase_class1 <- plyr::count(PPase_class, vars =c("Protein.ID","class"))%>%arrange(freq)
Count_PPase_class1$class <- factor(Count_PPase_class1$class, levels=c("PPM/PP2C","PPP","PTP","Asp-dep"))
Count_PPase_class1$Protein.ID<- factor(Count_PPase_class1$Protein.ID, levels=Count_PPase_class1%>%arrange(freq)%>%pull(1))
## Extended Data Fig2b ##
ggplot(data=Count_PPase_class1, aes(x=class, y=freq,fill=freq,group=Protein.ID))+geom_bar(color="white",stat="identity",show.legend =T,linewidth=0.1)+
  scale_fill_viridis(discrete = F,option = "G",direction = -1,begin =0.1,end =0.95)+
  geom_text(position =position_stack(vjust = 0.5),aes(label=freq),colour="white",size =3.2)+#coord_flip()+
  labs(y='P_sites number',x='')+scale_y_continuous(position = "left")+
  theme(axis.text=element_text(size=12,colour = "black"), 
        axis.title=element_text(size=14),
        axis.line = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_blank(), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#ggsave("../Code/Figure/Github/PPase count_change.png",width=4,height =4)
dev.off()

### legend of Extended Data Fig2b ##
legend <- ggplot(data=Count_PPase_class1, aes(x=class, y=freq,fill=Protein.ID,group=Protein.ID))+geom_bar(color="white",stat="identity",show.legend =T,linewidth=0.1)+
  scale_fill_viridis(discrete = T,option = "G",direction = -1,begin =0.1,end =0.95)+
  geom_text(position =position_stack(vjust = 0.5),aes(label=freq),colour="white",size =3.2)+#coord_flip()+
  labs(y='P_sites number',x='')+scale_y_continuous(position = "left")+
  theme(axis.text=element_text(size=12,colour = "black"), 
        axis.title=element_text(size=14),
        axis.line = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_blank(), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
deviat <- ggpubr::get_legend(legend)%>%as_ggplot()
#ggsave("../Code/Figure/Github/PPase count_change_legend.png",width=4,height =5)
### Then change color in Adobe illustator accordingly ####
dev.off()


## remove Solyc10g005640.3.1 by filter(WT_cluster!=syr1_cluster)
PPase_class <- rbind(PP2C[,c(1:2,5:6,109)],PPP[,c(1:2,5:6,109)],PTP[,c(1:2,5:6,109)],Asp[,c(1:2,5:6,109)])%>%filter(WT_cluster!=syr1_cluster)
Count_PPase_class2 <- plyr::count(PPase_class, vars =c("class","Protein.ID","WT_cluster"))%>%arrange(desc(Protein.ID))%>%filter(WT_cluster!="filtered")
ggplot(Count_PPase_class2,aes(WT_cluster,factor(Protein.ID,levels =unique(Protein.ID))))+facet_grid(~as.factor(class))+
  geom_point(aes(color=freq,size=freq),show.legend = F)+
  geom_text(aes(label=freq),colour="white",size =4)+scale_x_discrete(limits = c("1","2","3","4","5"))+
  scale_color_gradientn(colours=c("#b3eebe","#46bac2", "#371ea3"),guide=guide_colorbar(reverse=TRUE)) +
  scale_size_area(max_size =9) + 
  labs(y='',x='Phosphorylation pattern',color='P_sites number')+guides(size = "none")+
  theme(axis.title=element_text(size=14,colour = 'black'), 
        axis.text=element_text(size=14,colour = 'black'), 
        panel.background = element_rect(color='black',linewidth = 1,fill = "transparent"), 
        panel.grid.major.y = element_line(colour = "gray",linetype = "dotted"),
        panel.grid.major.x = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size=14,colour = 'black',hjust = 0.5, vjust = 0.5),
        legend.text = element_text(size=12,colour = 'black'),
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent",color = NA), # get rid of legend panel bg
        strip.background = element_rect(colour = "black",linewidth = 1),
        strip.text = element_text(colour = "gray25", face = "bold",size = 12),
        plot.background = element_rect(fill = "transparent", color = NA)) 
#ggsave("../Code/Figure/Github/PLL2 matters.png",width=8,height =4)
#PLL2 solyc06g076100
dev.off()

