## DEVELOP NEW PIPELINE FOR FACS ANALYSIS USING flowCore AND flowMeans

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

## IMPORT FCS FILE(S)
######################
filename <- "/Users/niclasthomas/Documents/flowR/FCS files/Treg exp_26_2014/Treg16102014_8009.fcs"
file <- read.FCS(filename, transformation=FALSE)

## TRANSFORM DATA FOR TWO-WAY PLOT
###################################
## Automatically estimate the logicle transformation based on the data
lgcl <- estimateLogicle(file, channels = c("FSC-A","SSC-A","FITC-A","PerCP-Cy5-5-A",
                                         "APC-A","APC-Cy7-A","Horizon v450-A",
                                         "Horizon v500-A","PE-A","PE-Texas Red-A","PE-Cy7-A"))

## Transform parameters using the estimated logicle transformation
trans.file <- transform(file, lgcl)

cells <- exprs(trans.file)[1:nrow(trans.file),1:11]
cells <- data.frame(cells)

xyplot( `Horizon v500-A` ~ `PE-A`, data=trans.file, smooth=FALSE, stat=T, xbin=128)

## BINNING RAW FCS DATA - COMPUTATIONALLY INEFFICIEN
#####################################################

breaks <- seq(0,4.5,by=1)
cat(paste("There will be",(length(breaks)+1)**dim(cells)[2],"bins"))

bins <- apply(cells, 2, cut, c(-Inf,breaks, Inf), labels=1:(length(breaks)+1))

appender <- function(x) paste(x, collapse='')

groups <- apply( bins, 1, appender )
table(groups)
max(table(groups))

## SETTING GATES FOR PLOTTING AND SUBSET DETERMINATION
######################################################
#n2gate <- norm2Filter("APC-Cy7-A","Horizon v500-A")
#c2gate <- curv2Filter("APC-Cy7-A","Horizon v500-A")
#c2gate <- curv2Filter("FSC-A","SSC-A")

#n2gate.results <- filter(trans.ex, n2gate)
#c2gate.results <- filter(ex,c2gate)

pca <- prcomp(~., data=cells, cor = TRUE, scale=T)
pca.df <- as.data.frame(pca$rotation)

pred <- predict(pca, cells)
pred <- as.data.frame(pred)

ggplot( pred, aes(PC1, PC2)) + geom_point()

## GATING STRATEGY USING flowMeans
##################################
fm.res <- flowMeans(trans.file, c("FSC-A","SSC-A","FITC-A",
                                "PerCP-Cy5-5-A",
                                "APC-A","APC-Cy7-A","Horizon v450-A",
                                "Horizon v500-A","PE-A","PE-Texas Red-A","PE-Cy7-A"
                                ),MaxN=10)

plot(trans.file[,c(1:11)], fm.res, c(forward,
                                   side,
                                   cd127,
                                   cd4,
                                   cd8,
                                   cd45ro,
                                   cd3,
                                   cd62l,
                                   pd1,
                                   hladr,
                                   cd25)
     )

library(ggplot2)

ggplot(dataset2, aes(x, y)) +
  stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
  scale_fill_gradientn(colours=colfunc(400)) + 
  xlim(c(15, 155)) + ylim(c(130, 270)) +
  geom_density2d(colour="black", bins=10) +
  geom_point() + 
  geom_text(aes(label=contesto), size=3, hjust=-.25, vjust=.75) +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10)) +
  theme(legend.title=element_blank())
