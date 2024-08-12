rm(list = ls())
setwd("~/Desktop/Dissertation/Publications/PLL5 paper/Manuscript/Source data/Code/")
getwd()
library(tidyverse)
library(ggplot2)
library(ggpubr)
#Solyc06g076100.3.1 S142-151-160 PLL2 
#The "X4.2_symbol" from phytozome ITAG3.2 was based on the Arabidopsis annotation(Best-hit-arabi-name), that's why it was called "PLL5" 

setup_combine <- readRDS("../Code/Data/setup_WTsyr1.rds")
df <- readRDS("../Code/Data/Summary_df.rds")

## P_sites measurement plot 
POI_df<- df%>%filter(Protein.ID=="Solyc06g076100.3.1")%>%select(c(1,21,37:108))%>%
  pivot_longer(cols = Intensity.WT_S01:Intensity.syr1_S456,names_to = "Samples", values_to = "Intensity")%>%
  mutate(genotype=word(word(Samples,1,sep = fixed("_")),2,sep = fixed(".")),
         site=as.numeric(word(word(ppep.site,2,sep =fixed(" ")),1,sep =fixed(";"))),
         STY_site=paste0(Amino.acid,word(word(ppep.site,2,sep =fixed(" ")),1,sep =fixed(";"))))%>%arrange(factor(site))
POI_df$minutes <- setup_combine[POI_df$Samples,"time.point"]
POI_df$genotype <- factor(POI_df$genotype ,levels=c("WT","syr1"))
POI_df$site <- factor(POI_df$site,levels=(sort(unique(POI_df$site))))
ggboxplot(POI_df, x = "genotype", y = "Intensity",color = "genotype", palette =c("#dc1400", "#00AFBB"),add = "jitter", shape = "genotype",
          facet.by = c("site"),panel.labs=list(site =unique(POI_df$STY_site)),bxp.errorbar=T,legend = "right",#scales = "free_y",
          title = unique(substr(POI_df$ppep.site,1,18)),subtitle =unique(df[POI_df$ppep.site,"X4.2_symbol"]),
          xlab = FALSE,ylab =c("Normalized intensity"))+theme(title = element_text(size=10,color = "black"))+
  geom_text(data=plyr::count(POI_df%>%filter(Intensity!="NA"), vars =c("site","genotype"))%>%mutate(Intensity=POI_df%>%filter(Intensity!="NA")%>%group_by(site,genotype)%>%summarise(low_bar = min(Intensity)-1)%>%pull(low_bar)),aes(label=paste0("N=",freq)),colour="black",parse=FALSE,size =3.5)


#### box plot:the whole set ######
ggboxplot(POI_df, x = "minutes", y = "Intensity",color = "genotype",add = ("jitter"),shape = "genotype",bxp.errorbar=T,
          panel.labs=list(site =unique(POI_df$STY_site)),facet.by = c("genotype","STY_site"),
          title = unique(substr(POI_df$ppep.site,1,18)),subtitle ="",
          xlab = c("Time (min)"),ylab ="Phosphorylation intensity",na.rm=T)+scale_color_manual(values=c("#dc1400","#00AFBB"))+guides(size=F,shape=F,color=F)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+
  geom_text(data=plyr::count(POI_df%>%filter(Intensity!="NA"), vars =c("genotype","minutes","STY_site"))%>%
  mutate(Intensity=POI_df%>%filter(Intensity!="NA")%>%group_by(genotype,minutes,STY_site)%>%
  summarise(low_bar = min(Intensity)-0.2)%>%pull(low_bar)),aes(label=paste0("N=",freq)),colour="black",parse=FALSE,size =3.5,nudge_x = 0)+
  theme(axis.text = element_text(size=10,color = "black"),title = element_text(size=10,color = "black"),
        plot.title = element_text(size=10,color = "black",hjust = 0),
        panel.background = element_rect(fill = "transparent"),strip.text = element_text(size=16,color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.line = element_line(colour = "black"),plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background =element_blank(),
        legend.background = element_rect(fill = "transparent",color = NA),legend.key = element_blank(),
        legend.box.background = element_rect(fill = "transparent",color = NA),legend.text=element_text(size=10))
#ggsave("../Code/Figure/Github/PLL2_3 sites.png",width=8,height=5)
dev.off()


## PLL2 S151 ##
int_combine <- readRDS("../Code/Data/int_norm.rds")%>%select(2:73)
ID <- "Solyc06g076100.3.1 151"
ID2 <- unique(df[POI_df$ppep.site,"X4.2_symbol"])
plot.data <- data.frame(minutes = as.numeric(as.character(setup_combine$time.IDs)), 
                        Intensity = as.numeric(int_combine[ID,]),
                        genotype = setup_combine[,"genotype"])#%>%filter(genotype=="WT")
## box plot: one by one
ggboxplot(plot.data, x = "minutes", y = "Intensity",color = "genotype",add = ("jitter"),shape = "genotype",bxp.errorbar=T,
          title = ID,facet.by = c("genotype"),#legend = "right",scales = "free_y",
          xlab = c("Time (min)"),ylab =F,na.rm=T)+#scale_y_continuous(limits=c(21, 25))+#scale_color_manual(values=c("#00AFBB"))+ 
  scale_color_manual(values=c("#dc1400","#00AFBB"))+scale_linetype_manual(values=c("dashed"))+
  stat_summary(aes(x = minutes, y = Intensity, color = genotype, group = genotype),fun=mean, geom="line",na.rm = TRUE,show.legend = F,linewidth =0.2,alpha=0.8)+
  geom_text(data=plyr::count(plot.data%>%filter(Intensity!="NA"), vars =c("genotype","minutes"))%>%
              mutate(Intensity=plot.data%>%filter(Intensity!="NA")%>%group_by(genotype,minutes)%>%summarise(low_bar = min(Intensity)-0.2)%>%pull(low_bar)),
            aes(label=paste0("N=",freq)),colour="black",parse=FALSE,size =3.5,nudge_x = 1)+
  theme(axis.text = element_text(size=10,color = "black"),title = element_text(size=10,color = "black"),
        plot.title = element_text(size=10,color = "black",hjust = 0),
        panel.background = element_rect(fill = "transparent"),strip.text = element_text(size=16,color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.line = element_line(colour = "black"),plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background =element_blank(),
        legend.background = element_rect(fill = "transparent",color = NA),
        legend.key = element_blank(),
        legend.box.background = element_rect(fill = "transparent",color = NA),
        legend.text=element_text(size=10))+guides(size=F,shape=F,color=F)
#ggsave("../Code/Figure/Github/S151_WTandsyr1.png",width=5,height =3)


