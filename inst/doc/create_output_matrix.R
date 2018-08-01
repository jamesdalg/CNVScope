## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = '~')
if(Sys.info()['sysname']=="Windows"){groupdir<-"W:/"} else {groupdir<-"/data/CCRBioinfo/"}

## ----echo=F,warning=F,message=F------------------------------------------
library(HiCNV)

## ------------------------------------------------------------------------
nbl_input_matrix<-readRDS("NBLTCGA_merged_df_aggregated_by_bin_fixed_comparisonv4.rds")
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

## ---- echo=F-------------------------------------------------------------
nbl_result_matrix<-readRDS("nbl_result_matrix_full.rds")

## ------------------------------------------------------------------------
nbl_result_matrix[1:5,1:5]

