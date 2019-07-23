## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

knitr::opts_knit$set(root.dir = '.')
library(CNVScope)

## ----eval=F,echo=T-------------------------------------------------------
#  library(CNVScope)

## ----eval=F,echo=T-------------------------------------------------------
#  if(!dir.exists("extracted_prad_data")){dir.create("extracted_prad_data")
#  untar("gdc_download_prad.tar.gz",exdir = "extracted_prad_data")}
#  tcga_files_prad<-list.files(path = "extracted_prad_data",pattern=glob2rx("*.tsv"),recursive=T,full.names = T)
#  print(tcga_files_prad)
#  

## ----eval=F,echo=T-------------------------------------------------------
#  sample_aggregated_segvals_output_full_prad<-formSampleMatrixFromRawGDCData(tcga_files = tcga_files_prad,format = "TCGA",binsize=1e6)
#  saveRDS(sample_aggregated_segvals_output_full_prad,"PRAD_sample_matched_input_matrix.rds")
#  

