## OPENCYTO PACKAGE PROCESSING FACS FILES
#########################################

rm(list=ls())

library(openCyto)
library(flowCore)

## READ IN T CELL GATING TEMPLATES
##################################
tcell.gating.template.control <- gatingTemplate("/home/niclas/Dropbox/facs/flowR/tcellgatingtemplate_control.csv") # control gating template
tcell.gating.template <- gatingTemplate("/home/niclas/Dropbox/facs/flowR/tcell_gatingtemplate.csv") # sample gating template


## READ IN FCS FILES, SPLITTING FILES INTO ISOTYPE CONTROL AND OTHER FILES
##########################################################################
path <- "/home/niclas/Documents/facstest/"
fcsFiles <- list.files(path)

iso.control <- grep("isotype",fcsFiles,ignore.case=TRUE,value=TRUE) # isotype control fcs file
samples <- grep("isotype",fcsFiles,ignore.case=TRUE,value=TRUE,invert=TRUE) # sample fcs files

samplesFlowSet <- read.flowSet(paste(path,samples,sep=""),alter.names=TRUE)
controlFlowSet <- read.flowSet(paste(path,iso.control,sep=""),alter.names=TRUE)
sampleNames(samplesFlowSet)<- c("baseline","week12")
sampleNames(controlFlowSet)<- "isotype.control"

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
gating(tcell.gating.template, samples.gating.set, mc.cores=2, parallel_type = "multicore")

## GET POPULATION STATS
#######################
summary.stats.control <- getPopStats(control.gating.set[[1]])
summary.stats.control
summary.stats.samples <- getPopStats(samples.gating.set[[1]])
summary.stats.samples

## USE PD1 GATES FROM ISOTYPE CONTROL ON SAMPLES
################################################

## GATE CD4 SUBSETS
for(i in list("CD4+CD8-/PD1","CD4+CD8-/CD25","CD4+CD8-/HLADR","CD4+CD8-/CD127")){
  try(Rm(i,samples.gating.set)) #remove the old gate
  gate <- getGate(control.gating.set, i) #get the new one
  add(samples.gating.set,gate$`isotype.control`, parent='CD4+CD8-') #add it
}
## GATE CD8 SUBSETS
for(i in list("CD4-CD8+/PD1","CD4-CD8+/CD25","CD4-CD8+/HLADR","CD4-CD8+/CD127")){
  try(Rm(i,samples.gating.set)) #remove the old gate
  gate <- getGate(control.gating.set, i) #get the new one
  add(samples.gating.set,gate$`isotype.control`, parent='CD4-CD8+') #add it
}
## GATE CD4 MEMORY SUBSETS
for(i in list(
  "CD62L+CD45RO+/PD1","CD62L+CD45RO+/CD25","CD62L+CD45RO+/HLADR","CD62L+CD45RO+/CD127",
  "CD62L+CD45RO-/PD1","CD62L+CD45RO-/CD25","CD62L+CD45RO-/HLADR","CD62L+CD45RO-/CD127",
  "CD62L-CD45RO+/PD1","CD62L-CD45RO+/CD25","CD62L-CD45RO+/HLADR","CD62L-CD45RO+/CD127",
  "CD62L-CD45RO-/PD1","CD62L-CD45RO-/CD25","CD62L-CD45RO-/HLADR","CD62L-CD45RO-/CD127")){
  try(Rm(i,samples.gating.set)) #remove the old gate
  gate <- getGate(control.gating.set, "CD4+CD8-/") #get the new one
  add(samples.gating.set,gate$`isotype.control`, parent=strsplit(i,"/")[[1]][1]) #add it
}
## GATE CD8 MEMORY SUBSETS
for(i in list(
  "CD62L+CD45RO+/PD1","CD62L+CD45RO+/CD25","CD62L+CD45RO+/HLADR","CD62L+CD45RO+/CD127",
  "CD62L+CD45RO-/PD1","CD62L+CD45RO-/CD25","CD62L+CD45RO-/HLADR","CD62L+CD45RO-/CD127",
  "CD62L-CD45RO+/PD1","CD62L-CD45RO+/CD25","CD62L-CD45RO+/HLADR","CD62L-CD45RO+/CD127",
  "CD62L-CD45RO-/PD1","CD62L-CD45RO-/CD25","CD62L-CD45RO-/HLADR","CD62L-CD45RO-/CD127")){
  try(Rm(i,samples.gating.set)) #remove the old gate
  gate <- getGate(control.gating.set, i) #get the new one
  add(samples.gating.set,gate$`isotype.control`, parent='CD62L+CD45RO-') #add it
}

recompute(samples.gating.set,alwaysLoadData=TRUE)
gating(tcell.gating.template,samples.gating.set)

## PLOT RESULTS
###############
#png('gating-strategy.png',width=5, height=5, units="in", res=500)
#plot(control.gating.set)
#dev.off()

plotGate(control.gating.set[[1]],"CD4+CD8-/PD1",smooth=FALSE,default.y="SSC.A",xbin=0)
plotGate(samples.gating.set[[1]],"CD4+CD8-/CD62L+CD45RO-/PD1",smooth=FALSE,default.y="SSC.A",xbin=0)

## GET DATA FOR UNDERLYING SUBSET
#################################
#tnaive <- getData(control.gating.set[[1]], "/nonDebris/Lymphocytes/CD3/CD4/PD1")
#tnaive<- as.data.frame(exprs(tnaive))
#tnaive <- tnaive[,-grep("Time",colnames(tnaive))]



i="CD62L+CD45RO-/PD1"
try(Rm(i,samples.gating.set)) #remove the old gate
gate <- getGate(control.gating.set, "CD4+CD8-/PD1")
add(samples.gating.set,gate$`isotype.control`, parent='CD62L+CD45RO-') #add it
