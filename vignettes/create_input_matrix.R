## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

knitr::opts_knit$set(root.dir = '.')
library(CNVScope)

## ----eval=F,echo=T-------------------------------------------------------
#  library(CNVScope)

## ----eval=F,echo=T-------------------------------------------------------
#  if(!dir.exists("extracted_nbl_data")){dir.create("extracted_nbl_data")}
#  untar("gdc_download_20180801_160142.tar.gz",exdir = "extracted_nbl_data")
#  tcga_files_nbl<-list.files(path = "extracted_nbl_data",pattern=glob2rx("*NormalVsPrimary.tsv"),recursive=T,full.names = T)
#  print(tcga_files_nbl)
#  

## ----eval=F,echo=T-------------------------------------------------------
#  sample_aggregated_segvals_output_full<-formSampleMatrixFromRawGDCData(tcga_files = tcga_files_nbl,format = "TARGET")
#  saveRDS(sample_aggregated_segvals_output_full,"NBL_sample_matched_input_matrix.rds")
#  

