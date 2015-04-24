## LOAD GATED FLOW SET & PROJECT MARKER FLUORESCENCE
## ON TO tSNE
####################################################
####################################################

rm(list=ls())

library(openCyto)
library(flowCore)
library(ggplot2)
library(Rtsne)
library(gridExtra)

setwd("/Users/niclasthomas/Dropbox/pid/bcellgatingsets")
gs_list <- lapply(list.files(),function(this_folder){
  load_gs(this_folder)
})

## DEFINE WHAT SUBSETS ARE DESIRED FOR ANALYSIS
getNodes(gs_list[[1]],order="bfs")
subsets <- getNodes(gs_list[[1]],order="bfs")[c(8)]

## GET ANTIBODY NAMES **bit of a hack - could improve this bit
metadata <- getData(gs_list[[1]])[[1]]
ab.names <- unlist(keyword(metadata,c("$P7S","$P8S","$P9S","$P10S","$P11S","$P12S","$P13S","$P14S","$P15S","$P16S","$P17S")))

## GET DATA FOR UNDERLYING SUBSET
#################################
get.underlying.data <- function(subsets,gating_list){
  total.data <- c()
  for (i in c(1:length(subsets))){
    flow.data <- lapply(gating_list,getData,y=subsets[i])
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
my.data <- c()
for (k in c(1:length(gs_list))){
  gs <- gs_list[[k]]
  nodes <- getNodes(gs,order="bfs")
  if (subsets %in% nodes){
    combined.data <- get.underlying.data(subsets,gs)
    combined.data <- cbind(rep(list.files()[k],dim(combined.data)[1]),combined.data)
    colnames(combined.data) <- c("sample",colnames(combined.data)[-1])
    my.data <- rbind(my.data,combined.data)
  }
}

## tSNE ON SUBSET POPULATION
############################
set.seed(6)
pred <- Rtsne(as.matrix(my.data[-c(1,2)]),theta=0.9)
pred <- as.data.frame(pred)

## 2D GGPLOT OF PC1 AND PC2
###########################
plot.data <- cbind(my.data[,c(1,2)],pred)

g <- ggplot( plot.data, aes(Y.1,Y.2,colour=sample)) +
  geom_point(size=1)+
  xlab("")+
  ylab("")+
  theme_bw()+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  guides(shape = guide_legend(override.aes = list(size=4)))
g
# ggsave(file=paste("/Users/niclasthomas/Dropbox/pid/figures/","compare_samples.png",sep=""),g,scale=1,width=10,height=10,dpi=300)
