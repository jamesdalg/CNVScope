# bin_data<-data.table::fread("ftp://helix.nih.gov/pub/dalgleishjl/plotly_dashboard_ext/binfile.txt")
# bin_data_t<-t(bin_data)
# colnames(bin_data_t)<-bin_data$probe
# bin_data_t<-
# bin_data_small1<-bin_data[1:100,]
# bin_data_small2=bin_data[1:100,]
# bin_data_small_t<-as.data.frame(t(bin_data[1:100,2:ncol(bin_data)]))
# colnames(bin_data_small_t)<-bin_data_small$probe
# #colnames()
# lm(data=bin_data_small_t,bin_data_small_t$chr10_118440805_118590868 ~ bin_data_small_t$chr10_118590869_118890996)
# #It's much easier to handle data with many rows than columns, hence, I'd go with having a long vs wide dataframe for this.
# calcRegForBins<-function(binmatrix)
# {
#   bin.pairs<-expand.grid(1:length(chromosomes),1:length(chromosomes))
# }
library(rslurm)
library(parallel)
calcVecLMs<-function(bin_data)
{
  if(dim(bin_data)[1]>dim(bin_data)[2]){bin_data<-t(bin_data)}
  bin_data_df<-as.data.frame(bin_data)
  bin.pairs<-expand.grid(1:ncol(bin_data),1:ncol(bin_data))
  neglogpvalues<-mcmapply(x=bin.pairs[,1],mc.cores=detectCores(),y=bin.pairs[,2],function(x,y) {-log(summary(lm(bin_data_df[,y]~bin_data_df[,x]))$coefficients[2,4])})
  output_matrix<-matrix(neglogpvalues,nrow=max(dim(bin_data)))
  return(output_matrix)
}
registerDoMC()
bin_data<-matrix(runif(189*14567),ncol=189)
start<-proc.time()
full_matrix<-calcVecLMs(bin_data[1:2000,])
proc.time()-start
#slurm version
library(rslurm)
bin_data<-matrix(runif(189*14567),ncol=189)
if(dim(bin_data)[1]>dim(bin_data)[2]){bin_data<-t(bin_data)}
bin_data_df<-as.data.frame(bin_data)
bin.pairs<-expand.grid(1:ncol(bin_data),1:ncol(bin_data))
pars<-
sjob <- slurm_apply(function(x,y) {-log(summary(lm(bin_data_df[,y]~bin_data_df[,x]))$coefficients[2,4])},
                    x=bin.pairs[,1],
y=bin.pairs[,2],
                    jobname = 'bin_data_linreg',
                    nodes = 56, cpus_per_node = 2, submit = FALSE)