## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
if(Sys.info()['sysname']=="Windows"){groupdir<-"W:/"} else {groupdir<-"/data/CCRBioinfo/"}

knitr::opts_knit$set(root.dir = paste0(groupdir,'dalgleishjl/hicnv/vignette/'))
setwd(paste0(groupdir,'dalgleishjl/hicnv/vignette/'))
library(HiCNV)

## ----eval=F--------------------------------------------------------------
#  library(HiCNV)

