## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

knitr::opts_knit$set(root.dir = '.')
library(magrittr)

## ----nbl_data_files_png,echo=F------------------------------------------------
knitr::include_graphics("NBL_data_files.png", dpi = 10)

## ----download_from_cart,echo=F------------------------------------------------
knitr::include_graphics("download_from_cart.png", dpi = 10)

## ----nbl_files,eval=F,echo=T--------------------------------------------------
#  if(!dir.exists("extracted_nbl_data")){dir.create("extracted_nbl_data")}
#  untar("gdc_download_20180801_160142.tar.gz",exdir = "extracted_nbl_data")
#  target_files_nbl<-list.files(path = "extracted_nbl_data",pattern=glob2rx("*NormalVsPrimary.tsv"),recursive=T,full.names = T)
#  print(target_files_nbl)
#  

## ----eval=F,echo=T------------------------------------------------------------
#  sample_aggregated_segvals_output_full<-CNVScope::formSampleMatrixFromRawGDCData(tcga_files = target_files_nbl,format = "TARGET")
#  saveRDS(sample_aggregated_segvals_output_full,"NBL_sample_matched_input_matrix.rds")
#  

## ---- eval=F,echo=T-----------------------------------------------------------
#  nbl_custom_input_matrix<-CNVScope::formSampleMatrixFromRawGDCData(tcga_files = target_files_nbl,
#  format = "custom",binsize = 1e6,freadskip = 14,parallel=F,debug=F,
#  sample_pat = "(?<=30-)(.*?)(?=_)",sample_col = "sample",chrlabel=">chr",
#  startlabel = "begin",endlabel = "end",cnlabel = "relativeCvg")
#  saveRDS(nbl_custom_input_matrix,"NBL_custom_sample_matched_input_matrix.rds")
#  

## ---- echo=T,eval=F-----------------------------------------------------------
#  nbl_custom_input_matrix_hd<-CNVScope::formSampleMatrixFromRawGDCData(tcga_files = target_files_nbl,
#  format = "custom",binsize = 2.5e5,freadskip = 14,parallel=T,debug=F,
#  sample_pat = "(?<=30-)(.*?)(?=_)",sample_col = "sample",chrlabel=">chr",
#  startlabel = "begin",endlabel = "end",cnlabel = "relativeCvg")
#  saveRDS(nbl_custom_input_matrix_hd,"NBL_custom_sample_matched_input_matrix_2.5e5binsize_parallel.rds")
#  

## ---- echo=T,eval=F-----------------------------------------------------------
#  nbl_custom_input_matrix_ld<-CNVScope::formSampleMatrixFromRawGDCData(tcga_files = target_files_nbl,
#  format = "custom",binsize = 1e8,freadskip = 14,parallel=F,debug=F,
#  sample_pat = "(?<=30-)(.*?)(?=_)",sample_col = "sample",chrlabel=">chr",
#  startlabel = "begin",endlabel = "end",cnlabel = "relativeCvg")
#  saveRDS(nbl_custom_input_matrix,"NBL_custom_sample_matched_input_matrix_1e8binsize.rds")
#  

