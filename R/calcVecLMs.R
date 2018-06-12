#'  Create a linear regression matrix.
#'
#' Creates a matrix of linear regression p-values, log transformed from every combination of columns in the parent matrix.
#' @keywords lm linear regression matrix
#' @import parallel rslurm
#' @param bin_data The parent matrix, with columns to have linear regression performed on them.
#' @param use_slurm Paralleize over a number of slurm HPC jobs? If false, the program will simply run locally.
#' @param slurmjob the slurm job object produced by rslurm::slurm_apply(), after running the function initially.
#' @return The output matrix, or if using slurm, the slurm job object (which should be saved as an rds file and reloaded when creating the output matrix).
#' 
#' @examples
#' #small example
#' #bin_data<-matrix(runif(5*5),ncol=5)
#' #full_matrix<-calcVecLMs(bin_data)
#' #real data example
#' #sample_aggregated_segvals<-readRDS("../NBLTCGA_merged_df_aggregated_by_bin_fixed_comparisonv4.rds")
#' @export
#' 
#' 
generateSignedNegLogPvalMat<-function(x,y)
{
  lmsum<-summary(lm(unlist(bin_data_df[,y])~unlist(bin_data_df[,x])));
  return(-log(lmsum$coefficients[2,4])*sign(lmsum$coefficients[2,1]))
  }  
calcVecLMs<-function(bin_data,use_slurm=F,job_finished=F,slurmjob=NULL,n_nodes=NULL,cpus_on_each_node=1,memory_per_node="2g",walltime="4:00:00",n_cores=1)
{
  #if(dim(bin_data)[1]<dim(bin_data)[2]){bin_data<-t(bin_data)}
  bin_data_df<-as.data.frame(bin_data)
  bin.pairs<-expand.grid(1:ncol(bin_data),1:ncol(bin_data))
  colnames(bin.pairs)<-c("x","y")
  if(!use_slurm){
    if(Sys.info()['sysname']!="Windows") {
      ncores=1
    } else {ncores=n_cores}
    if(Sys.info()['sysname']!="Windows")    { 
      neglogpvalues<-mcmapply(mc.cores = ncores,x=bin.pairs[,1],y=bin.pairs[,2],FUN = function(x,y) {lmsum<-summary(lm(unlist(bin_data_df[,y])~unlist(bin_data_df[,x]))); return(-log(lmsum$coefficients[2,4])*sign(lmsum$coefficients[2,1])) } ) } else {
      neglogpvalues<-mapply(x=bin.pairs[,1],y=bin.pairs[,2],FUN = function(x,y) {lmsum<-summary(lm(unlist(bin_data_df[,y])~unlist(bin_data_df[,x]))); return(-log(lmsum$coefficients[2,4])*sign(lmsum$coefficients[2,1])) } )
    }
    output_matrix<-matrix(neglogpvalues,ncol=ncol(bin_data))
    output_matrix[is.infinite(output_matrix)]<-max(output_matrix[!is.infinite(output_matrix)])
    return(output_matrix)
  }
  if(use_slurm){
    if(job_finished)
    {
      output_matrix<-matrix(get_slurm_out(slurmjob),ncol=ncol(bin_data))
      output_matrix[is.infinite(unlist(output_matrix))]<-max(unlist(output_matrix[!is.infinite(unlist(output_matrix))]))
      colnames(output_matrix)<-colnames(bin_data)
      rownames(output_matrix)<-colnames(bin_data)
      return(output_matrix)
    } else {
      if(is.null(n_nodes)){n_nodes<-ncol(bin_data_df)/2}
      lm_test_sjob <- slurm_apply(f = function(x,y) {lmsum<-summary(lm(unlist(bin_data_df[,y])~unlist(bin_data_df[,x]))); return(-log(lmsum$coefficients[2,4])*sign(lmsum$coefficients[2,1])) },params =  bin.pairs, jobname = 'test_apply',
                                  nodes = n_nodes, cpus_per_node = cpus_on_each_node, submit = T,slurm_options = list(partition="ccr,norm,quick",mem=memory_per_node,time=walltime,cpus_per_task=cpus_on_each_node))
      return(lm_test_sjob)
    }
  }
}