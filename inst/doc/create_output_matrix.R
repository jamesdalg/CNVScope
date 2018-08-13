## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = '~')
if(Sys.info()['sysname']=="Windows"){groupdir<-"W:/"} else {groupdir<-"/data/CCRBioinfo/"}
setwd("~")
#nbl_input_matrix<-readRDS("NBLTCGA_merged_df_aggregated_by_bin_fixed_comparisonv4.rds")
nbl_input_matrix<-readRDS("NBL_sample_matched_input_matrix.rds")

#nbl_result_matrix<-readRDS("nbl_result_matrix_full.rds")


## ----echo=F,warning=F,message=F------------------------------------------
library(HiCNV)

## ------------------------------------------------------------------------
nbl_input_matrix[1:5,1:5]

## ---- eval=F-------------------------------------------------------------
#  library(parallel)
#  nbl_slurm_object_test_zero_removed<-calcVecLMs(bin_data = as.data.frame(t(nbl_input_matrix))[,colSums(as.data.frame(t(nbl_input_matrix)))>0],use_slurm = T,n_nodes = 975,memory_per_node = "32g",walltime = "04:00:00",n_cores = 2,cpus_on_each_node = 2,job_finished = F,slurmjob = NULL)

## ---- eval=F-------------------------------------------------------------
#  
#  saveRDS(nbl_slurm_object_test_zero_removed,"nbl_slurm_object_test_zero_removed.rds")

## ---- eval=F-------------------------------------------------------------
#  nbl_result_matrix<-matrix(get_slurm_out(nbl_slurm_object_test_zero_removed),ncol=ncol(as.data.frame(t(nbl_input_matrix))[,colSums(as.data.frame(t(nbl_input_matrix)))>0]))
#  saveRDS(nbl_result_matrix,"nbl_result_matrix_full.rds")
#  

## ---- echo=F,include=T---------------------------------------------------
nbl_result_matrix<-readRDS("nbl_result_matrix_full.rds")
nbl_result_matrix_sign_corrected<-readRDS("nbl_result_matrix_sign_corrected.rds")

## ------------------------------------------------------------------------
nbl_result_matrix[1:5,1:5]

## ------------------------------------------------------------------------
nbl_result_matrix[1:5,1:5]

## ----eval=F--------------------------------------------------------------
#  nbl_result_matrix_sign_corrected<-postProcessLinRegMatrix(input_matrix = nbl_input_matrix,LM_mat = nbl_result_matrix,cor_type = "pearson",inf_replacement_val = 300)

## ------------------------------------------------------------------------
nbl_result_matrix_sign_corrected[1:5,1:5]
ComplexHeatmap::Heatmap(signedRescale(as.matrix(nbl_result_matrix_sign_corrected)),col = circlize::colorRamp2(c(0,0.5,1),c("blue","white","red")),cluster_rows = F,cluster_columns = F,show_heatmap_legend = F,show_column_names = F,show_row_names = F)

## ----eval=F--------------------------------------------------------------
#  if(!dir.exists("nbl_matrix_set")){dir.create("nbl_matrix_set")}
#  #setwd("nbl_matrix_set")
#  createChromosomalMatrixSet(whole_genome_mat=nbl_result_matrix_sign_corrected,output_dir="nbl_matrix_set",prefix="nbl_")

## ------------------------------------------------------------------------
setwd("nbl_matrix_set")
list.files()

