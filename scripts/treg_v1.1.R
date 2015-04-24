## OPENCYTO PACKAGE PROCESSING FACS FILES
#########################################

rm(list=ls())

library(openCyto)
library(flowCore)

setwd("/Users/niclasthomas/Dropbox/pid/")

healthy.control <- "/Users/niclasthomas/Documents/cvid-fcs-files/Exp_9_2014/Treg exp_9_2014/Treg_8029.fcs" 
fcsFile <- "/Users/niclasthomas/Documents/cvid-fcs-files/Exp_9_2014/Treg exp_9_2014/Treg_1039.fcs"

## READ IN T CELL GATING TEMPLATES
##################################
tcell.gating.template.control <- gatingTemplate(paste(getwd(),"/gatingtemplates/","gatingtemplatecontrol_treg.csv",sep="")) # control gating template
tcell.gating.template <- gatingTemplate(paste(getwd(),"/gatingtemplates/","gatingtemplate_treg.csv",sep="")) # sample gating template

## READ IN FCS FILES, SPLITTING FILES INTO ISOTYPE CONTROL AND OTHER FILES
##########################################################################
samplesFlowSet <- read.flowSet(fcsFile,alter.names=TRUE)
controlFlowSet <- read.flowSet(healthy.control,alter.names=TRUE)

sampleNames(controlFlowSet) <- "healthy.control"
sampleNames(samplesFlowSet) <- "1123"

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

## APPLY CD45RA GATE FROM CD4 ON TREG AND TCONV
###############################################
add(samples.gating.set, getGate(samples.gating.set, "CD4/CD45RA"), parent="Treg") #add it
add(samples.gating.set, getGate(samples.gating.set, "CD4/CD45RA"), parent="Tconv") #add it
try(Rm("CD4/CD45RA",samples.gating.set)) #remove the old gate

## USE GATES FROM ISOTYPE CONTROL ON SAMPLES
############################################
aliases <- list("Treg/CTLA4+","Treg/KI67+","Treg/ICOS+","Treg/CD39+","Treg/CD73+",
                "Treg/CTLA4-","Treg/KI67-","Treg/ICOS-","Treg/CD39-","Treg/CD73-")
parents <- list("Treg","Treg","Treg","Treg","Treg",
                "Treg","Treg","Treg","Treg","Treg")

for (i in c(1:length(aliases))){
  try(Rm(aliases[[i]],samples.gating.set)) #remove the old gate
  gate <- getGate(control.gating.set, aliases[[i]]) #get the new one
  add(samples.gating.set,gate$`healthy.control`, parent=parents[[i]]) #add it
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
#plot(samples.gating.set)
#dev.off()

# pdf('gating-strategy.pdf',width=5, height=5)
plotGate(samples.gating.set[[1]],"Treg",smooth=FALSE,xbin=0,default.y="SSC.A")
 # dev.off()
#plotGate(samples.gating.set, "Treg.PD1", type = "densityplot", stack = T)

## NORMALISE GATING SET
#######################
#my_target_sample <- sampleNames(samples.gating.set)[1]
#G.norm <- normalize(data=samples.gating.set)

## SAVE GATING SET TO DISK
##########################
save_gs(samples.gating.set,path=paste(getwd(),"/gatingsets/",sampleNames(samplesFlowSet),sep=""))

