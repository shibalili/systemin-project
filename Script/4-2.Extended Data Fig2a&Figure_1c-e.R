rm(list = ls())
setwd("~/Desktop/Dissertation/Publications/PLL5 paper/Manuscript/Source data/Code/")
getwd()
library(tidyverse)
library(stringr)
library(factoextra)
library(reshape2)# melt

## Optimal k ##
Smean_df <- readRDS("../Code/Data/kmc_ori_df.rds")
df <- t(scale(t(na.omit(Smean_df[,1:6]))))%>%as.data.frame()
set.seed(123)
fviz_nbclust(x = df,FUNcluster = kmeans,method = "gap_stat",diss =get_dist(df,method ="pearson"),nboot =100,verbose = TRUE,k.max =10)+labs(subtitle = "Gap Statistic Method")
#ggsave(filename = "../Code/Figure/Github/Gap Statistic Method.png",width =5,height =4)

## k=5 ##
Smean_df <- readRDS("../Code/Data/kmc_imp10_df.rds")
df <- t(scale(t(na.omit(Smean_df[,1:6]))))%>%as.data.frame()
final <- eclust(df,k=5, FUNcluster="kmeans",nstart =100,iter.max = 100, graph = FALSE,seed = 123)# kmeans/"euclidean"only
a<- data.frame(df,cluster= as.character(final$cluster),genotype=str_split_fixed(rownames(df),"_",2)[,2])
a$genotype <- factor(a$genotype, levels = c("WT","syr1"))
table(a$genotype)
a$cluster <- factor(a$cluster, levels = c("1","2","3","4","5"))
table(a$cluster,a$genotype)
a$ppep.site <- word(rownames(a),1,sep = fixed("_"))
colnames(a)[1:6] <- c("0","1","2","5","15","45")
#saveRDS("../Code/Data/k5_imp10.rds",object = a)

## Figure 1d ##
a <- readRDS("../Code/Data/k5_imp10.rds")
Order_a <- a[order(a$cluster,a$genotype),]
#saveRDS("../Code/Data/k5_imp10_Order_a.rds",object = Order_a)
distribution <- as.data.frame(Order_a%>%group_by(cluster=as.character(cluster))%>%dplyr::select(cluster,genotype,ppep.site)%>%pivot_wider(names_from=genotype, values_from=cluster))
distribution$Protein.ID <- substr(distribution$ppep.site,1,18)
rownames(distribution) <- distribution$ppep.site
distribution[,c("WT","syr1")][is.na(distribution[,c("WT","syr1")])] <- "filtered" 
#saveRDS("../Code/Data/k5_imp10_dist.rds",object = distribution)

rem <- distribution%>%filter(WT==3|WT==2,syr1==3|syr1==2)%>%pull(ppep.site) # remove syr1 cluster 2 or 3
Order_a <- Order_a%>%filter(!ppep.site%in%rem)
distribution <- as.data.frame(Order_a%>%group_by(cluster=as.character(cluster))%>%dplyr::select(cluster,genotype,ppep.site)%>%pivot_wider(names_from=genotype, values_from=cluster))
distribution$Protein.ID <- substr(distribution$ppep.site,1,18)
rownames(distribution) <- distribution$ppep.site
distribution[,c("WT","syr1")][is.na(distribution[,c("WT","syr1")])] <- "filtered" 
table(distribution$WT) # cluster 2:463 & cluster 3 573
#saveRDS("../Code/Data/k5_imp10_dist_no23.rds",object = distribution)

distribution <- readRDS("../Code/Data/k5_imp10_dist_no23.rds")%>%filter(WT!=syr1) # remove similar profiles shared in both WT and syr1
table(distribution$WT,distribution$syr1)
Order_a_rem <- Order_a[Order_a$ppep.site%in%rownames(distribution),]
int_norm <- readRDS("../Code/Data/int_norm.rds")
a <- Order_a_rem[Order_a_rem$ppep.site%in%rownames(int_norm),]%>%select(1:9)
b <- melt(a,id.vars = c("cluster","genotype","ppep.site"),variable.name = c("time_point"), value.name = "norm_intensity")
se <- function(x) sqrt(var(x)/length(x))#SE calculation
Percluster <- b %>% group_by(cluster, genotype, time_point) %>% summarize(mean=median(norm_intensity), SE=se(norm_intensity))
Percluster$cluster <- factor(Percluster$cluster, levels=c("1","2","3","4","5"))
ggplot(data = b, aes(x = time_point, y = norm_intensity,group=ppep.site))+
  geom_line(aes(color=genotype),linewidth =0.6,alpha=0.03,show.legend = F)+labs(x = "Time (min)",y= "Normalized intensity")+
  facet_grid(as.factor(genotype)~as.factor(cluster))+
  geom_line(Percluster,mapping=aes(y =mean,group=interaction(cluster,genotype)),alpha=0.5)+
  theme(axis.text = element_text(colour = "black",size=13),
        axis.title=element_text(size=16,colour = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA),
        strip.text= element_text(size =14,colour = "black"),
        panel.grid = element_blank(),strip.background = element_blank(),
        legend.position="top",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(fill = "transparent",color = NA),
        legend.title=element_text(size=22,face = "bold"),
        legend.text=element_text(size=22))+
  geom_text(data=plyr::count(a, vars =c("genotype","cluster")),aes(x=1.6, y=1.9, label=paste0("N=",freq)),colour="black", inherit.aes=FALSE, parse=FALSE,size = 4)
#ggsave("../Code/Figure/Github/k5_imp10.png",width =8,height =3.8)

## Figure 1e : change L48 with command below ##
distribution <- readRDS("../Code/Data/k5_imp10_dist_no23.rds")%>%filter(WT!=syr1,WT==3|WT==2) 
#ggsave("../Code/Figure/Github/k5_imp10_k=23.png",width =8,height =3.8)

dev.off()


## Figure 1c: stacking plot ##
distribution <- readRDS("../Code/Data/k5_imp10_dist_no23.rds")%>%filter(WT!=syr1)
Order_a <- readRDS("../Code/Data/k5_imp10_Order_a.rds")
Order_a_rem <- Order_a[Order_a$ppep.site%in%rownames(distribution),]
int_norm <- readRDS("../Code/Data/int_norm.rds")
a <- Order_a_rem[Order_a_rem$ppep.site%in%rownames(int_norm),]%>%select(1:9)
count_by_treat <- NULL
for (i in 1:5) {num <-data.frame(a%>%filter(cluster==i)%>%count(genotype),cluster=i)
count_by_treat <- rbind(count_by_treat,num)}
count_by_treat <- merge(count_by_treat,count_by_treat%>%group_by(cluster)%>%summarise(sum=sum(n)))
## y percentage 
ggplot(count_by_treat, aes(fill=genotype, y=n/sum, x=cluster)) + geom_bar(position="stack", stat="identity",show.legend = T,alpha=0.8)+
  scale_x_continuous(breaks=seq(1,7,1))+ylab("P_sites ratio")+ggtitle("Psites distribution in 5 clusters")+
  scale_fill_manual(name="Genotype",values= c("#dc1400",  "#00AFBB"),labels=c("WT","syr1"))+
  geom_text(aes(label=n),size=4,color="white",position = position_stack(vjust = 0.5))+#theme_minimal()+
  theme(plot.title = element_text(colour = "black",size=12),
        axis.title=element_text(size=12,colour = "black"),
        axis.text = element_text(colour = "black",size=10),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(fill = "transparent",color = NA),
        legend.title=element_text(size=12),
        #legend.position="top",
        legend.text=element_text(size=10))
ggsave("../Code/Figure/Github/stackplot_k5_imp10.png",width =3.5,height =3)
dev.off()







