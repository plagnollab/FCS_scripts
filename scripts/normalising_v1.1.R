###################################################
### code chunk number 1: loadGvHD
###################################################
library(flowStats)
data(ITN)

###################################################
### code chunk number 2: transform
###################################################
wf <- workFlow(ITN)
tl <- transformList(colnames(ITN)[3:7], asinh, transformationId="asinh")
add(wf, tl)

###################################################
### code chunk number 5: variation
###################################################
pars <- colnames(Data(wf[["base view"]]))[c(3,4,5,7)]
print(densityplot(PatientID~., Data(wf[["TCells+"]]), channels=pars, groups=GroupID,
                  scales=list(y=list(draw=F)), filter=lapply(pars, curv1Filter),
                  layout=c(4,1)))

###################################################
### code chunk number 6: norm
###################################################
norm <- normalization(normFun=function(x, parameters, ...)
  warpSet(x, parameters, ...),
  parameters=pars,
  arguments=list(grouping="GroupID", monwrd=TRUE),
  normalizationId="Warping")
add(wf, norm, parent="TCells+")


###################################################
### code chunk number 7: normPlot
###################################################
print(densityplot(PatientID~., Data(wf[["Warping"]]), channels=pars, groups=GroupID,
                  scales=list(y=list(draw=F)), filter=lapply(pars, curv1Filter),
                  layout=c(4,1)))