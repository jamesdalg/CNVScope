calcVecLMs<-function(bin_data,use_slurm=F,job_finished=F,slurmjob=NULL,n_nodes=NULL,cpus_on_each_node=1,memory_per_node="2g",walltime="4:00:00")
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
      output_matrix<-matrix(get_slurm_out(slurmjob),ncol=ncol(bin_data))
      output_matrix[is.infinite(unlist(output_matrix))]<-max(unlist(output_matrix[!is.infinite(unlist(output_matrix))]))
      colnames(output_matrix)<-colnames(bin_data)
      rownames(output_matrix)<-colnames(bin_data)
      return(output_matrix)
    } else {
      if(is.null(n_nodes)){n_nodes<-ncol(bin_data_df)/2}
      lm_test_sjob <- slurm_apply(function(x,y) -log(summary(lm(unlist(bin_data_df[,y])~unlist(bin_data_df[,x])))$coefficients[2,4]), bin.pairs, jobname = 'test_apply',
                                  nodes = n_nodes, cpus_per_node = cpus_on_each_node, submit = T,slurm_options = list(partition="ccr,norm,quick",mem=memory_per_node,time=walltime,cpus_per_task=cpus_on_each_node))
      return(lm_test_sjob)
    }
  }
}