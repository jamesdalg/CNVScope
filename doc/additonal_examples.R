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
chr_7_probdist <- CNVScope::calcCNVKernelProbDist(cor(chr_7_mat,use="pairwise.complete.obs"),parallel=F)$percentile_matrix
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


## ----census_data,eval=F--------------------------------------------------
#  census_data <- readRDS(system.file("plotly_dashboard_ext","censushg19.rds",package = "CNVScopePublicData"))
#  census_data[census_data@seqnames %in% "chr7"] %>% sort() %>% tibble::as_tibble() %>% janitor::clean_names() %>% dplyr::select(seqnames,start,end,gene_symbol,tumour_types_somatic,tumour_types_germline) %>% dplyr::filter(start>60e6,stringr::str_detect(string = tumour_types_somatic,pattern="AML") | stringr::str_detect(string = tumour_types_germline,pattern="AML"))

## ----blca_files,eval=F,echo=T--------------------------------------------
#  if(!dir.exists("extracted_blca_data")){dir.create("extracted_blca_data")
#  untar("gdc_download_blca.tar.gz",exdir = "./extracted_blca_data")}
#  tcga_files_blca<-list.files(path = "extracted_blca_data",pattern=glob2rx("*.tsv"),recursive=T,full.names = T)
#  print(tcga_files_blca)

## ----eval=F,echo=T-------------------------------------------------------
#  sample_aggregated_segvals_blca<-formSampleMatrixFromRawGDCData(tcga_files = tcga_files_blca,format = "TCGA",parallel=T)
#  saveRDS(sample_aggregated_segvals_blca,"blca_sample_matched_input_matrix.rds")

## ----blca_plots, eval=T,echo=T-------------------------------------------
sample_aggregated_segvals_blca<-readRDS("blca_sample_matched_input_matrix.rds")
invariant_bins<-which((sample_aggregated_segvals_blca[stringr::str_detect(rownames(sample_aggregated_segvals_blca),"chr17"),] %>% t() %>% as.data.frame() %>% sapply(sd))==0)
chr_17_mat<-sample_aggregated_segvals_blca[(stringr::str_detect(rownames(sample_aggregated_segvals_blca),"chr17") & rownames(sample_aggregated_segvals_blca) %in% setdiff(rownames(sample_aggregated_segvals_blca),names(invariant_bins))),] %>% t()

## ----chr17_cor-----------------------------------------------------------
chr_17_mat %>%  cor(use="pairwise.complete.obs",method="pearson") %>% 
  CNVScope::signedRescale(max_cap=1) %>%
  reshape2::melt()  %>%
  ggplot(aes(x=reshape2::colsplit(Var1,"_",c("chr","start","end"))$start,
             y=reshape2::colsplit(Var2,"_",c("chr","start","end"))$start,
             fill=value)) + geom_raster() +
  theme(axis.text.x = element_blank(),axis.text.y=element_blank(),axis.title = element_blank()) +
  ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.5, limits = c(0, 1))

## ----probdist_chr17------------------------------------------------------
chr_17_probdist <- CNVScope::calcCNVKernelProbDist(cor(chr_17_mat,use="pairwise.complete.obs"),parallel=F)$percentile_matrix
colnames(chr_17_probdist)<-colnames(chr_17_mat)
rownames(chr_17_probdist)<-colnames(chr_17_mat)
chr_17_js_breakpoints<-jointseg::jointSeg(chr_17_probdist,K=40)$bestBkp
chr_17_js_breakpoint_labels<-colnames(cor(chr_17_mat))[chr_17_js_breakpoints]
chr_17_js_breakpoint_labels


## ----breakpoint_plot_chr17,eval=F----------------------------------------
#  
#  breakpoint_plot_probdist <- chr_17_probdist %>% #  cor(use="pairwise.complete.obs",method="pearson") %>%
#      CNVScope::signedRescale(max_cap=1) %>%
#      reshape2::melt()  %>%
#    dplyr::mutate(col_pos=reshape2::colsplit(Var1,"_",c("chr","start","end"))$start,
#           row_pos=reshape2::colsplit(Var2,"_",c("chr","start","end"))$start,
#           rel_prob=value) %>%
#      ggplot(aes(x=col_pos,
#                 y=row_pos,
#                 fill=rel_prob)) + geom_raster() +
#      theme(axis.text.x = element_text(angle=90, hjust=1),axis.text.y=element_blank()) +
#      scale_x_continuous(breaks=reshape2::colsplit(chr_17_js_breakpoint_labels,"_",c("chr","start","end"))$start,labels=chr_17_js_breakpoint_labels) +
#      ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.5, limits = c(0, 1)) +
#    labs(x="col_pos",y="row_pos",value="Pearson Correlation:") + ggtitle("Chromosome 17 relationship probability") +
#    geom_contour(binwidth = .395, aes(z = value))
#  breakpoint_plot <- chr_17_mat %>%   cor(use="pairwise.complete.obs",method="pearson") %>%
#      CNVScope::signedRescale(max_cap=1) %>%
#      reshape2::melt()  %>%
#    dplyr::mutate(col_pos=reshape2::colsplit(Var1,"_",c("chr","start","end"))$start,
#           row_pos=reshape2::colsplit(Var2,"_",c("chr","start","end"))$start,
#           correlation=value) %>%
#      ggplot(aes(x=col_pos,
#                 y=row_pos,
#                 fill=correlation)) + geom_raster() +
#      theme(axis.text.x = element_text(angle=90, hjust=1),axis.text.y=element_blank()) +
#      scale_x_continuous(breaks=reshape2::colsplit(chr_17_js_breakpoint_labels,"_",c("chr","start","end"))$start,labels=chr_17_js_breakpoint_labels) +
#      ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.5, limits = c(0, 1)) +
#    labs(x="col_pos",y="row_pos",value="Pearson Correlation:") + ggtitle("Chromosome 17 linear relationship domains") +
#    geom_contour(binwidth = .395, aes(z = value))
#  breakpoint_plot_corr_diff <- ((chr_17_mat %>%   cor(use="pairwise.complete.obs",method="spearman"))-(chr_17_mat %>%   cor(use="pairwise.complete.obs",method="pearson"))) %>%
#      CNVScope::signedRescale(max_cap=1) %>%
#      reshape2::melt()  %>%
#    dplyr::mutate(col_pos=reshape2::colsplit(Var1,"_",c("chr","start","end"))$start,
#           row_pos=reshape2::colsplit(Var2,"_",c("chr","start","end"))$start,
#           corr_diff=value) %>%
#      ggplot(aes(x=col_pos,
#                 y=row_pos,
#                 fill=corr_diff)) + geom_raster() +
#      theme(axis.text.x = element_text(angle=90, hjust=1),axis.text.y=element_blank()) +
#      scale_x_continuous(breaks=reshape2::colsplit(chr_17_js_breakpoint_labels,"_",c("chr","start","end"))$start,labels=chr_17_js_breakpoint_labels) +
#      ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.5, limits = c(0, 1)) +
#    labs(x="col_pos",y="row_pos",value="Pearson Correlation:") + ggtitle("Chromosome 17 nonlinear (red) relationship regions, inferred by nonlinear-linear correlation difference") +
#    geom_contour(binwidth = .395, aes(z = value))
#  
#  breakpoint_plot
#  breakpoint_plot_probdist
#  breakpoint_plot_corr_diff
#  

## ----plotly_blca,eval=F--------------------------------------------------
#  library(plotly)
#  breakpoint_plot %>% plotly::ggplotly()

## ----3D_blca,eval=F------------------------------------------------------
#  chr_17_long <- chr_17_mat %>%   cor(use="pairwise.complete.obs",method="pearson") %>%
#      CNVScope::signedRescale(max_cap=1) %>%
#      reshape2::melt()  %>%
#    dplyr::mutate(col_pos=as.numeric(reshape2::colsplit(Var1,"_",c("chr","start","end"))$start),
#           row_pos=as.numeric(reshape2::colsplit(Var2,"_",c("chr","start","end"))$start),
#           correlation=value) %>% dplyr::select(col_pos,row_pos,correlation)
#  plot_ly(data = chr_17_long, x=chr_17_long$col_pos,y=chr_17_long$row_pos,z=chr_17_long$correlation,color=c(0,0.5,1),colors=colorRamp(c("blue","white","red")),intensity=chr_17_long$correlation,type = "mesh3d")

## ----skcm_files,eval=F,echo=T--------------------------------------------
#  if(!dir.exists("extracted_skcm_data")){dir.create("extracted_skcm_data")}
#  untar("gdc_download_skcm.tar.gz",exdir = "./extracted_skcm_data")
#  tcga_files_skcm<-list.files(path = "extracted_skcm_data",pattern=glob2rx("*.tsv"),recursive=T,full.names = T)
#  print(tcga_files_skcm)

## ----eval=F,echo=T-------------------------------------------------------
#  #ptm <- proc.time()
#  #doMC::registerDoMC()
#  #doParallel::registerDoParallel()
#  sample_aggregated_segvals_skcm<-formSampleMatrixFromRawGDCData(tcga_files = tcga_files_skcm,format = "TCGA",parallel = T)
#  #proc.time() - ptm
#  saveRDS(sample_aggregated_segvals_skcm,"skcm_sample_matched_input_matrix.rds")

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

