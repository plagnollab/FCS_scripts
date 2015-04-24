## LOAD GATED FLOW SET & PROJECT MARKER FLUORESCENCE
## ON TO PRINCIPAL COMPONENTS
####################################################
####################################################

rm(list=ls())

library(openCyto)
library(flowCore)
library(ggplot2)
library(gridExtra)
library(Rtsne)

setwd("/Users/niclasthomas/Dropbox/pid/bcellgatingsets")
# gs_list <- lapply(list.files(),function(this_folder){
#   load_gs(this_folder)
# })
filename <- "1230"
gs <- load_gs(filename)

## DEFINE WHAT SUBSETS ARE DESIRED FOR ANALYSIS
getNodes(gs,order="bfs")
subsets <- getNodes(gs,order="bfs",path="auto")[c(5:13)]

## GET ANTIBODY NAMES **bit of a hack - could improve this bit
metadata <- getData(gs)
ab.names <- unlist(keyword(metadata,c("$P7S","$P8S","$P9S","$P10S","$P11S","$P12S","$P13S","$P14S","$P15S","$P16S","$P17S")))

getPopStats(gs,"count")

## GET DATA FOR UNDERLYING SUBSET
#################################
get.underlying.data <- function(subsets,gs){
  total.data <- c()
  for (i in c(1:length(subsets))){
    flow.data <- getData(gs,subsets[i])
    subset.data <- data.frame()
    for (j in c(1:length(flow.data))){
      x <- flow.data[[j]]
      x <- x[,-grep("Time",colnames(x))]
      x <- x[,-grep("FSC",colnames(x))]
      x <- x[,-grep("SSC",colnames(x))]
      subset.data <- rbind(subset.data,exprs(x))
    }
    subset.data <- cbind(subsets[i],subset.data)
    colnames(subset.data) <- c("subset",colnames(subset.data)[-1])
    total.data <- rbind(total.data,subset.data)
  }
  return(total.data)
}

## LOOP OVER ALL GATING SETS IN FOLDER
######################################
#gs <- gs_list[[k]]
combined.data <- get.underlying.data(subsets,gs)

## tSNE ON SUBSET POPULATION
############################
set.seed(6)
#  size <- c()
#   for (i in subsets){
#     size <- c(size,length(which(combined.data[,c(1)]==i)))
#   }
#   sample.size <- min(n,min(size))
# my.sampled.data <- c()
# for (j in subsets){
#   chosen.samples <- sample(which(combined.data[,c(1)]==j),sample.size,replace=T)
#   my.sampled.data <- rbind(my.sampled.data,combined.data[chosen.samples,])
# }
pred <- Rtsne(as.matrix(combined.data[,-c(1)]),theta=0.5,check_duplicates=FALSE)
pred <- as.data.frame(pred)

## 2D GGPLOT OF PC1 AND PC2
###########################
pred <- cbind(combined.data[,c(1)],pred)
colnames(pred) <- c("subset",colnames(pred)[-1])
g <- ggplot( pred, aes(x=Y.1,y=Y.2,colour=subset)) +
  geom_point(aes(colour=subset),size=1)+
  #theme_bw()+
  labs(colour="Subset")+
  xlab("")+
  ylab("")+
  guides(colour = guide_legend(override.aes = list(size=4)))
  #scale_colour_brewer(palette="Set1")
print(g)
ggsave(file=paste("/Users/niclasthomas/Dropbox/pid/bcellfigures/tsne-nosampling/",filename,".png",sep=""),g,scale=1,width=10,height=10,dpi=300)

