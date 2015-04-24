## PERFORM PCA ON SUBSET FREQUENCIES (+/- MFI)
## FROM CVID PATIENTS
##############################################
rm(list=ls())

library(ggplot2)

setwd("/Users/niclasthomas/Dropbox/pid/")
input <- "data/B cell summary.csv"

## READ DATA
############
data <- read.csv(input,sep=",",header=T)
data <- data[complete.cases(data),]
data <- data[-grep("frozen",data$Sample),]

patient.names <- sapply(strsplit(as.character(data$Sample),"_"), function(x) x[[2]])
patient.names <- sapply(strsplit(as.character(patient.names),".",fixed=TRUE), function(x) x[[1]])

## COLOUR POINTS BASED ON CLINICAL STATUS
#########################################
# colrejection <- c()
# for (i in c(1:length(truth))){
#   if (truth[i] == "R"){colrejection <- c(colrejection,2)}
#   if (truth[i] == "NR"){colrejection <- c(colrejection,19)}
# }

## PERFORM PCA
##############
pca <- prcomp(~., data=data[,-c(1)], cor = TRUE, scale=T)
pca.scores <- as.data.frame(pca$x)
pca.scores <- cbind(patient.names,pca.scores)

## PLOT PCA
###########
g <- ggplot(data=pca.scores,aes(x=PC1,y=PC2,label=patient.names))+
  geom_point(size=2)+
  geom_text(aes(label=patient.names),hjust=0,vjust=0,size=4)+
  xlab("PC1")+
  ylab("PC2")+
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(-4,4))+
  theme_bw()
g
ggsave(file=paste(getwd(),"/figures/","pca-Bcell.png",sep=""),g,scale=1)

## DETERMINE LOADINGS
#####################
# W1 <- pca$rotation[,1]
# w1.summary <- round(W1[order(W1)],digits=3)
# 
# W2 <- pca$rotation[,2]
# w2.summary <- round(W2[order(W2)],digits=3)
# 
# topn <- 5
# for (i in c(1:topn)){
#   print(w1.summary[i])
# }
# for (i in c(length(w1.summary):(length(w1.summary)-topn+1))){
#   print(w1.summary[i])
# }
# 
# topn <- 5
# for (i in c(1:topn)){
#   print(w2.summary[i])
# }
# for (i in c(length(w2.summary):(length(w2.summary)-topn+1))){
#   print(w2.summary[i])
# }
