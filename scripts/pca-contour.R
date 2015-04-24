## DEVELOP NEW PIPELINE FOR FACS ANALYSIS USING flowCore AND flowMeans
######################################################################

rm(list=ls())

## REFERENCE
#   FITC-A            CD127
#   PerCP-Cy5-5-A     HLA-DR
#   APC-A             CD62L
#   APC-Cy7-A         CD25
#   Horizon v450-A    CD8
#   Horizon v500-A    CD4
#   PE-A              PD-1
#   PE-Texas Red-A    CD45RO
#   PE-Cy7-A          CD3

## IMPORT LIBRARIES
####################
library(flowCore)
library(flowViz)
library(flowMeans)
library(ggplot2)
library(gtools)

## IMPORT FCS FILE(S)
#####################
myfiles <- list("/home/niclas/Documents/facsfiles/20130319/5048TCells20130319/T Cell Panel_Baseline G221711002504_32.fcs",
                "/home/niclas/Documents/facsfiles/20140519/5064 T Cells 20140519/T Cells_Baseline G221712000880_32.fcs"
)

allcells <- c()

for (i in myfiles){
  currentfile <- i
  file <- read.FCS(currentfile, transformation=FALSE)
  
  ## TRANSFORM DATA FOR TWO-WAY PLOT
  ###################################
  ## Automatically estimate the logicle transformation based on the data
  lgcl <- estimateLogicle(file, channels = c("FSC-A","SSC-A","FITC-A","PerCP-Cy5-5-A",
                                             "APC-A","APC-Cy7-A","Horizon v450-A",
                                             "Horizon v500-A","PE-A","PE-Texas Red-A","PE-Cy7-A"))
  
  ## Transform parameters using the estimated logicle transformation
  trans.file <- transform(file, lgcl)
  
  cells <- exprs(trans.file)[1:nrow(trans.file),1:11]
  
  allcells <- rbind(allcells,cells)
}
  
  ## DETERMINE PCA SCORES
  #######################
  pca <- prcomp(~., data=cells, cor = TRUE, scale=T)
  pca.df <- as.data.frame(pca$rotation)
  
  pred <- predict(pca, cells)
  pred <- as.data.frame(pred)
  
  # ggplot( pred, aes(PC1, PC2)) +
  #   geom_point() +
  #   xlim(c(-5, 5)) + ylim(c(-5, 5))
  
  #colfunc <- colorRampPalette(c("darkblue", "lightblue", "green", "yellow", "red"))
  
  dims <- 2
  dataset <- pred[,1:dims] # Select PCs 1, 2 and 3 only
  
  ## PLOT CONTOURS OF PCA
  #######################
  # ggplot(dataset, aes(PC1, PC2)) +
  #   stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
  #   scale_fill_gradientn(colours=colfunc(400)) +
  #   xlim(c(-5, 5)) + ylim(c(-5, 5)) +
  #   geom_density2d(colour="black", bins=5)
  
  breaks <- seq(-5,5,by=2)
  #cat(paste("There will be",(length(breaks)+1)**dim(dataset)[2],"bins"))
  
  bins <- apply(dataset, 2, cut, c(-Inf,breaks, Inf), labels=letters[1:(length(breaks)+1)])
  
  appender <- function(x) paste(x, collapse='')
  
  groups <- apply( bins, 1, appender )
  counts <- as.data.frame(table(groups),stringsAsFactors=FALSE)
  
  letset <-  letters[1:(length(breaks)+1)]
  perms <- permutations(n=length(letset),r=dims,v=letset,repeats.allowed=T)
  allgroups <- apply(perms,1,appender)
  toadd <- setdiff(allgroups,counts$groups)
  for (i in 1:length(toadd)){counts <- rbind(counts,c(toadd[i],0))}
  counts <- counts[order(counts$groups),]
  
  write.table( t(counts$Freq), file = paste(outpath,outfile,sep=""), append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
