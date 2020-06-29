#'  Create a linear regression matrix.
#'
#' Creates a matrix of linear regression p-values, log transformed from every combination of columns in the parent matrix.
#' @keywords lm linear regression matrix
#' @import parallel
#' @importFrom rslurm slurm_apply get_slurm_out
#' @param bin_data The parent matrix, with columns to have linear regression performed on them.
#' @param use_slurm Paralleize over a number of slurm HPC jobs? If false, the program will simply run locally.
#' @param slurmjob the slurm job object produced by rslurm::slurm_apply(), after running the function initially.
#' @param job_finished Are all the slurm jobs finished and the results need retrieving?
#' @param n_nodes the number of nodes used in your slurm job.
#' @param cpus_on_each_node The number of cpus used on each node
#' @param memory_per_node the amount of ram per node (e.g. "32g" or "2g")
#' @param walltime Time for job to be completed for SLURM scheduler in hh:mm:ss format. Defaults to 4h.
#' @return The output matrix, or if using slurm, the slurm job object (which should be saved as an rds file and reloaded when creating the output matrix).
#' @examples
#'
#' #small example
#' #bin_data<-matrix(runif(5*5),ncol=5)
#' foreach::registerDoSEQ()
#' #full_matrix<-suppressWarnings(calcVecLMs(bin_data))
#' #Please note that lm() will make a warning when there are two vectors that are too close 
#' #numerically (this will always happen along the diagonal).
#' #This is normal behavior and is controlled & accounted for using this function as well as
#' #the postProcessLinRegMatrix function (which converts the infinite values to a maximum).
#'
#' @export
calcVecLMs<-function(bin_data,use_slurm=F,job_finished=F,slurmjob=NULL,n_nodes=NULL,cpus_on_each_node=2,memory_per_node="2g",walltime="4:00:00")
{
  #if(dim(bin_data)[1]<dim(bin_data)[2]){bin_data<-t(bin_data)}
  bin_data_df<-as.data.frame(bin_data)
  bin.pairs<-expand.grid(1:ncol(bin_data),1:ncol(bin_data))
  colnames(bin.pairs)<-c("x","y")
  if(!use_slurm){
    neglogpvalues<-mcmapply(x=bin.pairs[,1],y=bin.pairs[,2],function(x,y) -log(summary(lm(unlist(bin_data_df[,y])~unlist(bin_data_df[,x])))$coefficients[2,4]) )
    output_matrix<-matrix(neglogpvalues,ncol=ncol(bin_data))
    output_matrix[is.infinite(output_matrix)]<-max(output_matrix[!is.infinite(output_matrix)])
    return(output_matrix)
  }
  if(use_slurm){
    if(job_finished)
    {
      output_matrix<-matrix(rslurm::get_slurm_out(slurmjob),ncol=ncol(bin_data))
      output_matrix[is.infinite(unlist(output_matrix))]<-max(unlist(output_matrix[!is.infinite(unlist(output_matrix))]))
      colnames(output_matrix)<-colnames(bin_data)
      rownames(output_matrix)<-colnames(bin_data)
      return(output_matrix)
    } else {
      if(is.null(n_nodes)){n_nodes<-ncol(bin_data_df)/2}
      lm_test_sjob <- rslurm::slurm_apply(function(x,y) -log(summary(lm(unlist(bin_data_df[,y])~unlist(bin_data_df[,x])))$coefficients[2,4]), bin.pairs, jobname = 'test_apply',
                                  nodes = n_nodes, cpus_per_node = cpus_on_each_node, submit = T,slurm_options = list(partition="ccr,norm,quick",mem=memory_per_node,time=walltime,cpus_per_task=cpus_on_each_node))
      return(lm_test_sjob)
    }
  }
}