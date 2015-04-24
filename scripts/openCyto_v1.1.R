## OPENCYTO PACKAGE PROCESSING FACS FILES
#########################################

rm(list=ls())

library(openCyto)
library(flowCore)

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

## READ IN T CELL GATING TEMPLATES
##################################
tcell.gating.template.control <- gatingTemplate("/home/niclas/Dropbox/facs/flowR/tcellgatingtemplate_control.csv") # control gating template
tcell.gating.template <- gatingTemplate("/home/niclas/Dropbox/facs/flowR/tcell_gatingtemplate.csv") # sample gating template

## COMPENSATION
###############
samplesFlowSet.comp <- fsApply(samplesFlowSet,function(frame){
  comp <- keyword(frame)$`SPILL`
  new_frame <- compensate(frame,comp)
  new_frame
})

controlFlowSet.comp <- fsApply(controlFlowSet,function(frame){
  comp <- keyword(frame)$`SPILL`
  new_frame <- compensate(frame,comp)
  new_frame
})

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
summary.stats <- getPopStats(control.gating.set[[1]])
summary.stats

## GET DATA FOR UNDERLYING SUBSET
#################################
#tnaive <- getData(control.gating.set[[1]], "/nonDebris/Lymphocytes/CD3/CD4/PD1")
#tnaive<- as.data.frame(exprs(tnaive))
#tnaive <- tnaive[,-grep("Time",colnames(tnaive))]


## GATE SAMPLES BASED ON ISOTYPE CONTROLS
#########################################
system.time(gating(tcell.gating.template, samples.gating.set, mc.cores = 4, parallel_type = "multicore", stop.at='CD4'))

## USE PD1 GATES FROM ISOTYPE CONTROL ON SAMPLES
################################################
Rm("PD1",samples.gating.set)
for(i in list("PD1")){
  try(Rm(i,samples.gating.set)) #remove the old gate
  gate <- getGate(control.gating.set, i) #get the new one
  add(samples.gating.set,gate$`isotype.control`, parent='CD4') #add it
}
recompute(samples.gating.set,alwaysLoadData=TRUE)
system.time(gating(tcell.gating.template,samples.gating.set))

## PLOT RESULTS
###############
png('gating-strategy.png',width=5, height=5, units="in", res=500)
plot(control.gating.set)
dev.off()

png('gated-samples-baseline.png',width=5, height=5, units="in", res=500)
plotGate(samples.gating.set[[1]],smooth=FALSE,xbin=0)
dev.off()

png('gated-samples-week12.png',width=5, height=5, units="in", res=500)
plotGate(samples.gating.set[[2]],smooth=FALSE,xbin=0)
dev.off()

png('gated-controls.png',width=5, height=5, units="in", res=500)
plotGate(control.gating.set[[1]],smooth=FALSE,xbin=0)
dev.off()

#plotGate(control.gating.set[[1]],"PD1",type="densityplot")
