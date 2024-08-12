rm(list = ls())
setwd("~/Desktop/Dissertation/Publications/PLL5 paper/Manuscript/Source data/Code/")
getwd()
library(tidyverse)
library(factoextra) #eclust #fviz_dend
int_norm <- readRDS("../Code/Data/int_norm.rds")%>%select(c(2:73))
setup <- readRDS("../Code/Data/setup_WTsyr1.rds")
rep_sum <- NULL
for (i in 1:6) {
  setup.new <- setup[setup[,5]=="WT"&setup[,3]==i,]
  s <- int_norm[,rownames(setup.new)]
  rep_sum<-cbind(rep_sum,apply(s, 1, function(x){return(sum(x, na.rm = TRUE))}))
}
colnames(rep_sum) <- c("WT_rep1","WT_rep2","WT_rep3","WT_rep4","WT_rep5","WT_rep6")
WT_sum_rep6 <- rep_sum

rep_sum <- NULL
for (i in 1:6) {
  setup.new <- setup[setup[,5]=="syr1"&setup[,3]==i,]
  s <- int_norm[,rownames(setup.new)]
  rep_sum<-cbind(rep_sum,apply(s, 1, function(x){return(sum(x, na.rm = TRUE))}))
}
colnames(rep_sum) <- c("syr1_rep1","syr1_rep2","syr1_rep3","syr1_rep4","syr1_rep5","syr1_rep6")
syr1_sum_rep6 <- rep_sum
df <- as.data.frame(cbind(WT_sum_rep6,syr1_sum_rep6))
#saveRDS("../Code/Data/Aggreg_WTsyr1.rds",object = df)

df <- readRDS("../Code/Data/Aggreg_WTsyr1.rds")
df1 <- t(na.omit(t(scale(t(df)))))
# Hierarchical clustering
hc.res <- eclust(df1, "hclust", k=2, hc_metric = "pearson", hc_method = "complete",graph = FALSE,seed=123)
# Visualize dendrograms
fviz_dend(hc.res,horiz = F,k_colors = c("#FC4E07","#00AFBB"),color_labels_by_k = TRUE,rect=T,rect_fill = T,lwd = 0.5,cex=0.8, main = "Dendrogram_experiment_setup")
ggsave("../Code/Figure/Github/Dendrogram_experiment_setup.png",width =5,height=4)
dev.off()
