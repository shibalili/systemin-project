rm(list = ls())
setwd("~/Desktop/Dissertation/Publications/PLL5 paper/Manuscript/Source data/Code/")
getwd()
library(hypeR)
library(tidyverse)
int_norm <- readRDS("../Code/Data/int_norm.rds")%>%mutate(Protein.ID=substr(ppep.site,1,18),.after = ppep.site)
back.ID <- unique(int_norm$Protein.ID) ##2155
distribution <- readRDS("../Code/Data/k5_imp10_dist.rds")
summary_df <- readRDS("../Code/Data/Summary_df.rds")
gl <- list()
for(x in 1:5){gl[[x]]<- distribution%>%filter(WT!=syr1,WT==x)%>%pull(Protein.ID)%>%unique()}
names(gl) <- c("k1","k2","k3","k4","k5")
lengths(gl)
anno <- readRDS("../Code/Data/ORA/BackID2155_Annotation_detail_bin1.rds")#Bin level:1 √ 
#anno <- readRDS("../Code/Data/ORA/BackID2155_Annotation_detail_bin1.1.rds")Bin level:1.1
IDs.sum <- apply(anno[,-c(1:4)], 2, sum)
str(IDs.sum)
relevant.bin <-names(IDs.sum)[IDs.sum>0]#how many bins(functional groups) in background
gl_bin <- list()
for(k in 1:5){
  gene.list <- gl[[k]]
  gene.set <- list() 
  anno.gene.list <- list() 
  for(i in 1:length(relevant.bin)){
    bin <- relevant.bin[i]
    gene.set[[i]] <- rownames(anno)[anno[,bin]==1]#as long as background is set, universe geneset is fixed
    anno.gene.list[[i]] <- intersect(gene.set[[i]], gene.list)}# anno genelist is varied due to diff clusters
  names(gene.set) <- relevant.bin #combine protein ID (from background list) with bincode
  names(anno.gene.list) <- relevant.bin
  gl_bin [[k]]<- names(anno.gene.list[lapply(anno.gene.list, length) > 0 ])}
names(gl_bin) <- names(gl)

b <- list() #map gene ID from background to binname
annotation <- readRDS("../Code/Data/ORA/Annotation source/ITAG3.2 mapman_Ensambl44.rds")
for(j in 1:length(gl_bin)){
  b[[j]] <- gene.set[gl_bin[[j]]]
  bin.name <- annotation[,2:3] %>% filter(BINCODE %in% gl_bin[[j]]) %>%unique()
  names(b[[j]]) <-paste0("bin ",bin.name[,1]," ",bin.name[,2])}
names(b) <- names(gl)

hyp_obj <- list()
for (m in 1:5) {hyp_obj[[m]] <- hypeR(signature=gl[[m]], genesets=b[[m]],test = "hypergeometric",background=length(back.ID))}
names(hyp_obj) <- names(gl) 
names(hyp_obj)

hyp_dots(hyp_obj[[2]],val = "pval",abrv = 50,title = "cluster2 bin1")
hyp_show(hyp_obj[[1]])#show table
dev.off()

###### organize df for better read or export #########
hyp_to_excel <- NULL
for (z in 1:5) {
  #hyp_to_excel$cluster=z
  df <- data.frame(hyp_obj[[z]]$data,cluster=z)
  hyp_to_excel <- rbind(hyp_to_excel,df)}
rownames(hyp_to_excel) <- NULL
y <- NULL
for (n in 1:5) {
  s <- strsplit(hyp_obj[[n]]$data[["hits"]], split = ",")
  y <-rbind(y,data.frame(sig_cluster=rep(names(hyp_obj)[n], sum(lengths(s))),
                         bin.name= rep(hyp_obj[[n]]$data[["label"]], sapply(s, length)),
                         pvalue= rep(hyp_obj[[n]]$data[["pval"]], sapply(s, length)),
                         Protein.ID = unlist(s)%>%gsub("\\s+","", .))) }# remove whitespace in a string
ppep.y1 <- NULL
for (r in 1:5) {
  site_kn <- summary_df%>%filter(WT_cluster!=syr1_cluster,WT_cluster==r)
  y_kn <- y%>%filter(sig_cluster==names(gl)[r])
  ppep.y <- subset(site_kn,site_kn$Protein.ID %in% y_kn$Protein.ID) ##map protein ID to ppep ID
  ppep.y1 <- rbind(ppep.y1,merge(y_kn,ppep.y,by ="Protein.ID"))%>%arrange(sig_cluster,pvalue)}
#saveRDS("../Code/Data/ORA/imp10_allBin1_k5-HA.rds",object=ppep.y1)
#saveRDS("../Code/Data/ORA/imp10_allBin1.1_k5-HA.rds",object=ppep.y1)



