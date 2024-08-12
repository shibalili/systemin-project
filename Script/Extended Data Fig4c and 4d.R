rm(list = ls())
setwd("~/Desktop/Dissertation/Publications/PLL5 paper/Manuscript/Source data/Code/")
getwd()
library(tidyverse)
library(ggplot2)
#Solyc07g017780.3.1 951 LHA4
#Solyc03g113405.1.1 780;949 LHA1 in ITAG3.2 which is no longer present in ITAG4 or 5  (Solyc03g113400.4.1 in ITAG4)
# The Phospho-site position needs to double check in the full length protein
# Often MaxQuant matches a site(especially conserved sites) to multiple phospho-peptide, the outcome of site position is not accurate in this case.
# In the case of penultimate(second to last) Thr of LHA1 and LHA4
# In "Flanking" column where the phospho-site is centered by 7aa(before) and 7aa(after)
# Solyc07g017780.3.1 951 LHA4 ETIQQHYTV______ 
# Solyc03g113405.1.1 780;949  ETIQQSYTV______


setup_combine <- readRDS("../Code/Data/setup_WTsyr1.rds")
df <- readRDS("../Code/Data/Summary_df.rds")
POI_df<- df%>%filter(Protein.ID=="Solyc03g113405.1.1")%>%select(c(1,21,37:108))%>%
  pivot_longer(cols = Intensity.WT_S01:Intensity.syr1_S456,names_to = "Samples", values_to = "Intensity")%>%
  mutate(genotype=word(word(Samples,1,sep = fixed("_")),2,sep = fixed(".")),
         site=as.numeric(word(word(ppep.site,2,sep =fixed(" ")),1,sep =fixed(";"))),
         STY_site=paste0(Amino.acid,word(word(ppep.site,2,sep =fixed(" ")),1,sep =fixed(";"))))%>%arrange(factor(site))
POI_df$minutes <- setup_combine[POI_df$Samples,"time.point"]
POI_df$genotype <- factor(POI_df$genotype ,levels=c("WT","syr1"))
POI_df$site <- factor(POI_df$site,levels=(sort(unique(POI_df$site))))

int_combine <- readRDS("../Code/Data/int_norm.rds")%>%select(2:73)
ID <- "Solyc03g113405.1.1 780;949"
ID2 <- ""
plot.data <- data.frame(minutes = as.numeric(as.character(setup_combine$time.IDs)), 
                        Intensity = as.numeric(int_combine[ID,]),
                        genotype = setup_combine[,"genotype"])#%>%filter(genotype=="WT")
ggplot(plot.data,aes(x = minutes, y = Intensity, color = genotype, group = genotype,linetype = genotype))+
  geom_point(na.rm = TRUE)+facet_wrap(vars(genotype))+stat_summary(fun=mean, geom="line",na.rm = TRUE)+
  scale_x_continuous(breaks=0:5, labels = c("0","1","2","5","15","45"))+#scale_y_continuous(limits = c(22,27))+
  labs(title =ID,subtitle = ID2)+ylab(label = "Phosphorylation intensity")+xlab("Time (min)")+
  scale_linetype_manual(values=c("solid", "dashed"))+scale_color_manual(values=c("#dc1400","#00AFBB"))+
  theme(axis.title = element_text(size=12,color = "black"),axis.line = element_line(colour = "black"),
        axis.text = element_text(size=8,color = "black"),
        plot.title = element_text(size=12,color = "black",hjust = 0),
        plot.subtitle = element_text(color = "black",size=10,hjust = 0),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"),panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_blank(),
        legend.background = element_rect(fill = "transparent",color = NA),legend.key = element_blank(),
        legend.box.background = element_rect(fill = "transparent",color = NA),legend.text=element_text(size=10))

#ggsave("../Code/Figure/Github/LHA4_17780 T951 box.png",width=5,height =3)
#ggsave("../Code/Figure/Github/LHA1_113400 T955 box.png",width=5,height =3)


dev.off()

