## ----setup, include=F----------------------------------------------------
knitr::opts_chunk$set(echo = F)
if(Sys.info()['sysname']=="Windows"){groupdir<-"W:/"} else {groupdir<-"/data/CCRBioinfo/"}

#knitr::opts_knit$set(root.dir = paste0(groupdir,"dalgleishjl/hicnv/hicnv/vignettes/")) #
#setwd("~")
#nbl_input_matrix<-readRDS("NBLTCGA_merged_df_aggregated_by_bin_fixed_comparisonv4.rds")
#getwd()
nbl_input_matrix<-readRDS("NBL_sample_matched_input_matrix.rds")

#nbl_result_matrix<-readRDS("nbl_result_matrix_full.rds")


## ----echo=T,warning=F,message=F------------------------------------------
library(CNVScope)

## ----echo=T--------------------------------------------------------------
nbl_input_matrix[1:5,1:5]

## ---- eval=F,echo=T------------------------------------------------------
#  library(parallel)
#  nbl_slurm_object_test_zero_removed<-calcVecLMs(bin_data =as.data.frame(t(nbl_input_matrix[which(rowSds(as.matrix(nbl_input_matrix))!=0.0),])),use_slurm = T,n_nodes = 975,memory_per_node = "32g",walltime = "04:00:00",cpus_on_each_node = 2,job_finished = F,slurmjob = NULL)
#  

## ---- eval=F,echo=T------------------------------------------------------
#  
#  saveRDS(nbl_slurm_object_test_zero_removed,"nbl_slurm_object_test_zero_removed.rds")

## ---- eval=F,echo=T------------------------------------------------------
#  library(matrixStats)
#  nbl_result_matrix<-matrix(as.numeric(unlist( get_slurm_out(nbl_slurm_object_test_zero_removed))),ncol=ncol(as.data.frame(t(nbl_input_matrix[which(rowSds(as.matrix(nbl_input_matrix))!=0.0),])) ) )
#  
#  saveRDS(nbl_result_matrix,"nbl_result_matrix_full.rds")
#  saveRDS(nbl_result_matrix[1:25,1:25],"nbl_result_matrix_small.rds")
#  

## ---- echo=F,include=T,eval=F--------------------------------------------
#  #source download locations
#  #download.file("https://github.com/jamesdalg/CNScope_public_data/blob/master/nbl_result_matrix_full.rds?raw=true","nbl_result_matrix_full.rds")
#  #download.file("https://github.com/jamesdalg/CNScope_public_data/blob/master/nbl_result_matrix_sign_corrected.rds","nbl_result_matrix_sign_corrected.rds")
#  #download.file("https://www.dropbox.com/s/sevuhos976c6guu/nbl_result_matrix_full.rds?dl=1","nbl_result_matrix_full.rds")
#  #download.file("https://www.dropbox.com/s/85hii5cd5epuuby/nbl_result_matrix_sign_corrected.rds?dl=1","nbl_result_matrix_sign_corrected.rds")
#  nbl_result_matrix<-readRDS("nbl_result_matrix_full.rds")
#  nbl_result_matrix_sign_corrected<-readRDS("nbl_result_matrix_sign_corrected.rds")

## ------------------------------------------------------------------------

## ----eval=T,echo=T-------------------------------------------------------
nbl_result_matrix_small<-readRDS("nbl_result_matrix_small.rds")
nbl_result_matrix_small[1:5,1:5]
nbl_result_matrix_sign_corrected<-postProcessLinRegMatrix(input_matrix = nbl_input_matrix[1:25,1:25],LM_mat = nbl_result_matrix_small,cor_type = "pearson",inf_replacement_val = 300)
nbl_result_matrix_sign_corrected[1:5,1:5]

## ----echo=T--------------------------------------------------------------
nbl_result_matrix_sign_corrected[1:5,1:5]
ComplexHeatmap::Heatmap(signedRescale(as.matrix(nbl_result_matrix_sign_corrected)),col = circlize::colorRamp2(c(0,0.5,1),c("blue","white","red")),cluster_rows = F,cluster_columns = F,show_heatmap_legend = F,show_column_names = F,show_row_names = F)

## ----eval=F,echo=T-------------------------------------------------------
#  if(!dir.exists("nbl_matrix_set")){dir.create("nbl_matrix_set")}
#  #setwd("nbl_matrix_set")
#  doMC::registerDoMC()
#  #use ONLY the whole matrix with chromosomes 1-X, not the small subset provided for documentation purposes.
#  createChromosomalMatrixSet(whole_genome_mat=nbl_result_matrix_sign_corrected,output_dir="nbl_matrix_set",prefix="nbl_")

## ----echo=T,eval=F-------------------------------------------------------
#  list.files("nbl_matrix_set")

