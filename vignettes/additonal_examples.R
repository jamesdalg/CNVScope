## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

knitr::opts_knit$set(root.dir = '.')
library(CNVScope)
options(scipen=999)
library(magrittr)

## ----aml_files,eval=F,echo=T---------------------------------------------
#  if(!dir.exists("extracted_aml_data")){dir.create("extracted_aml_data")}
#  untar("gdc_download_aml.tar.gz",exdir = "./extracted_aml_data")
#  target_files_aml<-list.files(path = "extracted_aml_data",pattern=glob2rx("*NormalVsPrimary.tsv"),recursive=T,full.names = T)
#  print(target_files_aml)

## ----eval=F,echo=T-------------------------------------------------------
#  sample_aggregated_segvals_aml<-formSampleMatrixFromRawGDCData(tcga_files = target_files_aml,format = "TARGET")
#  saveRDS(sample_aggregated_segvals_aml,"aml_sample_matched_input_matrix.rds")

## ----aml_plots, eval=T,echo=T--------------------------------------------
sample_aggregated_segvals_aml<-readRDS("aml_sample_matched_input_matrix.rds")
invariant_bins<-which((sample_aggregated_segvals_aml[stringr::str_detect(rownames(sample_aggregated_segvals_aml),"chr7"),] %>% t() %>% as.data.frame() %>% sapply(sd))==0)
chr_7_mat<-sample_aggregated_segvals_aml[(stringr::str_detect(rownames(sample_aggregated_segvals_aml),"chr7") & rownames(sample_aggregated_segvals_aml) %in% setdiff(rownames(sample_aggregated_segvals_aml),names(invariant_bins))),] %>% t()

## ----chr7_cor------------------------------------------------------------
chr_7_mat %>%  cor(use="pairwise.complete.obs",method="pearson") %>% 
  CNVScope::signedRescale(max_cap=1) %>%
  reshape2::melt()  %>%
  ggplot(aes(x=reshape2::colsplit(Var1,"_",c("chr","start","end"))$start,
             y=reshape2::colsplit(Var2,"_",c("chr","start","end"))$start,
             fill=value)) + geom_raster() +
  theme(axis.text.x = element_blank(),axis.text.y=element_blank(),axis.title = element_blank()) +
  ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.5, limits = c(0, 1))

## ----breakpoints---------------------------------------------------------
colnames(chr_7_mat)[CNVScope::getAsymmetricBlockIndices(cor(chr_7_mat,use="pairwise.complete.obs"))]
breakpoints<-colnames(chr_7_mat)[CNVScope::getAsymmetricBlockIndices(cor(chr_7_mat,use="pairwise.complete.obs"))] %>% stringr::str_split_fixed(string = .,pattern="_",n=3) %>% as.matrix() %>% .[,2] %>% as.numeric()
breakpoint_labels <- colnames(chr_7_mat)[CNVScope::getAsymmetricBlockIndices(cor(chr_7_mat,use="pairwise.complete.obs"))]
breakpoint_labels

## ----breakpoint_plot-----------------------------------------------------
chr_7_mat %>%  cor(use="pairwise.complete.obs",method="pearson") %>% 
    CNVScope::signedRescale(max_cap=1) %>%
    reshape2::melt()  %>%
    ggplot(aes(x=reshape2::colsplit(Var1,"_",c("chr","start","end"))$start,
               y=reshape2::colsplit(Var2,"_",c("chr","start","end"))$start,
               fill=value)) + geom_raster() +
    theme(axis.text.x = element_text(angle=90, hjust=1),axis.text.y=element_blank(),axis.title = element_blank()) +
    scale_x_continuous(breaks=breakpoints,labels=breakpoint_labels) +
    ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.5, limits = c(0, 1))


## ----probdist------------------------------------------------------------
chr_7_probdist <- CNVScope::calcCNVKernelProbDist(cor(chr_7_mat,use="pairwise.complete.obs"))$percentile_matrix
js_breakpoints<-jointseg::jointSeg(chr_7_probdist,K=20)$bestBkp
js_breakpoint_labels<-colnames(chr_7_mat)[js_breakpoints]


## ----plot_probdist-------------------------------------------------------
chr_7_probdist %>%  
#  CNVScope::signedRescale(max_cap=1) %>%
  reshape2::melt()  %>%
  ggplot(aes(x=Var1,
             y=Var2,
             fill=value)) + geom_tile() +
#  theme(axis.title = element_blank()) + #axis.text.x = element_blank(),axis.text.y=element_blank(),
    theme(axis.text.x = element_text(angle=90, hjust=1),
          axis.text.y = element_text(angle=0, hjust=1)
          ,axis.title = element_blank()) +
#    scale_x_continuous(breaks=js_breakpoints,labels=js_breakpoint_labels) +
#      scale_y_continuous(breaks=js_breakpoints,labels=js_breakpoint_labels) +
      scale_x_continuous(breaks=js_breakpoints,labels=js_breakpoint_labels) +
      scale_y_continuous(breaks=js_breakpoints,labels=js_breakpoint_labels) +

  ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.5, limits = c(0, 1))


