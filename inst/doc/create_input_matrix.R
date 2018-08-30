## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
if(Sys.info()['sysname']=="Windows"){groupdir<-"W:/"} else {groupdir<-"/data/CCRBioinfo/"}

knitr::opts_knit$set(root.dir = paste0(groupdir,'dalgleishjl/hicnv/vignette/'))
#setwd(paste0(groupdir,'dalgleishjl/hicnv/vignette/'))
library(HiCNV)

## ----eval=F--------------------------------------------------------------
#  library(HiCNV)

## ------------------------------------------------------------------------
if(!dir.exists("extracted_nbl_data")){dir.create("extracted_nbl_data")}
untar("gdc_download_20180801_160142.tar.gz",exdir = "extracted_nbl_data")
tcga_files_nbl<-list.files(path = "extracted_nbl_data",pattern=glob2rx("*NormalVsPrimary.tsv"),recursive=T,full.names = T)
print(tcga_files_nbl)


## ----eval=F--------------------------------------------------------------
#  sample_aggregated_segvals_output_full<-formSampleMatrixFromRawGDCData(tcga_files = tcga_files_nbl,format = "TARGET")
#  saveRDS(sample_aggregated_segvals_output_full,"NBL_sample_matched_input_matrix.rds")
#  

