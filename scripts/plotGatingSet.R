## LOAD GATED FLOW SET & PROJECT MARKER FLUORESCENCE
## ON TO PRINCIPAL COMPONENTS
####################################################
####################################################

library(openCyto)

setwd("/Users/niclasthomas/Dropbox/pid/treggatingsets")
path <- "/Users/niclasthomas/Dropbox/pid/tregfigures/gatedplots/"

gs_list <- lapply(list.files(),function(this_folder){
  load_gs(this_folder)
})

for (i in c(1:length(gs_list))){
  pdf(paste(path,'gatedPlots-',sampleNames(gs_list[[i]]),'.pdf',sep=""))
  plotGate(gs_list[[i]][[1]],smooth=FALSE,xbin=32,default.y="SSC.A",path="auto")
  dev.off()
}