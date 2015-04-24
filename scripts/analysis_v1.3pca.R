## LOAD GATED FLOW SET & PROJECT MARKER FLUORESCENCE
## ON TO PRINCIPAL COMPONENTS
####################################################
####################################################

rm(list=ls())

library(openCyto)
library(flowCore)
library(ggplot2)
library(rgl)
library(gridExtra)

setwd("/Users/niclasthomas/Dropbox/pid/testfolder")
gs_list <- lapply(list.files(),function(this_folder){
  load_gs(this_folder)
})

## DEFINE WHAT SUBSETS ARE DESIRED FOR ANALYSIS
getNodes(gs_list[[1]],order="bfs")
subsets <- getNodes(gs_list[[1]],order="bfs",path="auto")[c(4)]

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

for (k in c(1:length(gs_list))){
  
  gs <- gs_list[[k]]
  combined.data <- get.underlying.data(subsets,gs)
  
  ## PCA ON SUBSET POPULATION
  ###########################
  pca <- prcomp(~., data=combined.data[,-c(1)], cor = TRUE, scale=T)
  pred <- predict(pca, combined.data[,-c(1)])
  pred <- cbind(combined.data[,1],pred)
  pred <- as.data.frame(pred)
  colnames(pred) <- c("subset",colnames(pred)[-1])
  pred$subset <- as.factor(pred$subset)
  
  ## 2D GGPLOT OF PC1 AND PC2
  ###########################
  g <- ggplot(pred, aes(PC1, PC2)) +
    geom_point(aes(color=subset),size=0.2)+
    theme_bw()+
    scale_color_discrete(name="T Cell Subset",
                         breaks=c(1:length(subsets)),
                         labels=subsets)+
    guides(colour = guide_legend(override.aes = list(size=4)))
  print(g)
  #ggsave(file=paste("/Users/niclasthomas/Dropbox/pid/figures/",list.files()[k],".png",sep=""),g,scale=1,width=10,height=10,dpi=300)
  
  scores <- as.data.frame(pca$rotation)
  scores <- cbind(ab.names,scores)
  print(scores)
  
  ## SHOW PC SCORES
  #################
  p1 <- ggplot(data=scores,aes(x=ab.names,y=PC1))+
    geom_bar(stat="identity",colour="black")+
    xlab("")+
    ylim(c(-1,1))+
    theme(axis.text.x=element_text(angle=90, vjust=0.5, size=13,colour="black",hjust=1),
          axis.text.y=element_text(size=13,angle=90,colour="black"))
  p2 <- ggplot(data=scores,aes(x=ab.names,y=PC2))+
    geom_bar(stat="identity",colour="black")+
    xlab("")+
    ylim(c(-1,1))+
    theme(axis.text.x=element_text(angle=90, vjust=0.5, size=13,colour="black",hjust=1),
          axis.text.y=element_text(size=13,angle=90,colour="black"))
  #pdf(paste("/Users/niclasthomas/Dropbox/pid/figures/",list.files()[k],"-PCscores",".pdf",sep=""))
  #grid.arrange(p1,p2,ncol=2)
  #dev.off()
}
