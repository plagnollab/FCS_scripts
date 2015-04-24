## OPENCYTO PACKAGE PROCESSING FACS FILES
#########################################

rm(list=ls())

library(openCyto)
library(flowCore)
library(devtools)
dev_mode()
install_github("RGLab/flowStats")
library(flowStats)

setwd("/home/niclas/Dropbox/facs/flowR")

id <- "5002"
iso.control <- "/home/niclas/Documents/facsfiles/20140211/5002 T Cells 20140211/5002 Isotypes_D0 G221711001111_33.fcs" 
fcsFile <- "/home/niclas/Documents/facsfiles/20140211/5002 T Cells 20140211/T Cell Panel_D0 G221711001111_33.fcs"

## READ IN T CELL GATING TEMPLATES
##################################
tcell.gating.template.control <- gatingTemplate("tcellgatingtemplate_control.csv") # control gating template
tcell.gating.template <- gatingTemplate("tcellgatingtemplate.csv") # sample gating template

## READ IN FCS FILES, SPLITTING FILES INTO ISOTYPE CONTROL AND OTHER FILES
##########################################################################
samplesFlowSet <- read.flowSet(fcsFile,alter.names=TRUE)
controlFlowSet <- read.flowSet(iso.control,alter.names=TRUE)

sampleNames(controlFlowSet)<- "isotype.control"
sampleNames(samplesFlowSet)<- c(paste("patient",id,sep=""))

## COMPENSATION
###############
apply.compensation <- function(frame){
colnames(keyword(frame)$`SPILL`) <- gsub("-",".",colnames(keyword(frame)$`SPILL`))
colnames(keyword(frame)$`SPILL`) <- gsub(" ",".",colnames(keyword(frame)$`SPILL`))
comp <- keyword(frame)$`SPILL`
new_frame <- compensate(frame,comp)
new_frame
}

samplesFlowSet.comp <- fsApply(samplesFlowSet,apply.compensation)
controlFlowSet.comp <- fsApply(controlFlowSet,apply.compensation)

## REMOVE VARIABLES FSC, SSC .. FROM LOGICLE TRANSFORM
######################################################
vars <- colnames(samplesFlowSet.comp)
vars <- vars[-grep("FSC",vars)]
vars <- vars[-grep("SSC",vars)]
vars <- vars[-grep("Time",vars)]

## TRANSFORM DATA USING LOGICLE TRANSFORM
#########################################
lgcl <- estimateLogicle(samplesFlowSet.comp[[1]], channels = vars) # estimate parameters to logicle transform samples and control
samplesFlowSet.trans <- transform(samplesFlowSet.comp, lgcl)
controlFlowSet.trans <- transform(controlFlowSet.comp, lgcl)

## PERFORM AUTOMATED GATING
###########################
control.gating.set <- GatingSet(controlFlowSet.trans)
gating(tcell.gating.template.control, control.gating.set, mc.cores=2, parallel_type = "multicore")

samples.gating.set <- GatingSet(samplesFlowSet.trans)
gating(tcell.gating.template, samples.gating.set, mc.cores=2, parallel_type = "multicore", stop.at="Treg" )

## USE GATES FROM ISOTYPE CONTROL TO FIND TREG POPULATION
#########################################################
aliases <- list("CD4.CD25","CD4.CD127")
parents <- list("CD4+CD8-","CD4+CD8-")

for (i in c(1:length(aliases))){
  try(Rm(aliases[[i]],samples.gating.set)) #remove the old gate
  gate <- getGate(control.gating.set, aliases[[i]]) #get the new one
  add(samples.gating.set,gate$`isotype.control`, parent=parents[[i]]) #add it
}

recompute(samples.gating.set,alwaysLoadData=TRUE)
gating(tcell.gating.template,samples.gating.set)

## USE GATES FROM ISOTYPE CONTROL ON SAMPLES
############################################
aliases <- list("CD4.CD25","CD4.CD127","CD4.PD1","CD4.HLADR",
                "TN4.PD1","TEM4.PD1","TCM4.PD1","TES4.PD1",
                "TN4.CD25","TEM4.CD25","TCM4.CD25","TES4.CD25",
                "TN4.HLADR","TEM4.HLADR","TCM4.HLADR","TES4.HLADR",
                "TN4.CD127","TEM4.CD127","TCM4.CD127","TES4.CD127",
                
                "CD8.CD25","CD8.CD127","CD8.PD1","CD8.HLADR",
                "TN8.PD1","TEM8.PD1","TCM8.PD1","TES8.PD1",
                "TN8.CD25","TEM8.CD25","TCM8.CD25","TES8.CD25",
                "TN8.HLADR","TEM8.HLADR","TCM8.HLADR","TES8.HLADR",
                "TN8.CD127","TEM8.CD127","TCM8.CD127","TES8.CD127",
                
                "Treg.PD1","Treg.HLADR",
                "TregN.PD1","TregEM.PD1","TregCM.PD1","TregES.PD1",
                "TregN.HLADR","TregEM.HLADR","TregCM.HLADR","TregES.HLADR")

parents <- list("CD4+CD8-","CD4+CD8-","CD4+CD8-","CD4+CD8-",
                "TN4","TEM4","TCM4","TES4",
                "TN4","TEM4","TCM4","TES4",
                "TN4","TEM4","TCM4","TES4",
                "TN4","TEM4","TCM4","TES4",
                
                "CD4-CD8+","CD4-CD8+","CD4-CD8+","CD4-CD8+",
                "TN8","TEM8","TCM8","TES8",
                "TN8","TEM8","TCM8","TES8",
                "TN8","TEM8","TCM8","TES8",
                "TN8","TEM8","TCM8","TES8",
                
                "Treg","Treg",
                "TregN","TregEM","TregCM","TregES",
                "TregN","TregEM","TregCM","TregES")

for (i in c(1:length(aliases))){
  try(Rm(aliases[[i]],samples.gating.set)) #remove the old gate
  gate <- getGate(control.gating.set, aliases[[i]]) #get the new one
  add(samples.gating.set,gate$`isotype.control`, parent=parents[[i]]) #add it
}

recompute(samples.gating.set,alwaysLoadData=TRUE)
gating(tcell.gating.template,samples.gating.set)

## GET POPULATION STATS
#######################
summary.stats.control <- getPopStats(control.gating.set[[1]])
summary.stats.control
summary.stats.samples <- getPopStats(samples.gating.set[[1]])
summary.stats.samples

## PLOT RESULTS
###############
#png('gating-strategy.png',width=5, height=5, units="in", res=500)
#plot(control.gating.set)
#dev.off()

#plot(samples.gating.set)
#plotGate(control.gating.set[[1]],"CD8.CD62L",smooth=FALSE,default.y="SSC.A",xbin=0,type="densityplot")
#plotGate(samples.gating.set[[1]],"TN8",smooth=FALSE,default.y="SSC.A",xbin=0)

#plotGate(samples.gating.set, "Treg.PD1", type = "densityplot", stack = T)

## NORMALISE GATING SET
#######################
#my_target_sample <- sampleNames(samples.gating.set)[1]
#G.norm <- normalize(data=samples.gating.set)

## SAVE GATING SET TO DISK
##########################
save_gs(samples.gating.set,path=paste("/home/niclas/Documents/facs-results/patient",id,sep=""))
