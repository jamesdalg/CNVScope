#' Calculate the probability distribution of CNV concordance events with a fast kernel
#'
#' This function produces several matrices, including a Z-score matrix from a matrix of the same size and a percentile matrix of these Z-scores
#' @param submatrix A matrix of CNV data in an intrachromosomal region (e.g. chr1 vs chr1 or chr5 vs chr5)
#' @param win a window size for the matrix that calculates the windowed average using the kernel function
#' @keywords CNV kernel probability distribution concordance fast
#' @import ComplexHeatmap foreach doMC spatialfil 
#' @export
#' @examples
#' set.seed(303)
#' mat<-matrix(data=runif(n = 25),nrow=5,ncol=5)
#' calcCNVKernelProbDist(matrix)
#' mat_prob_dist<-calcCNVKernelProbDist(mat)
#' mat_prob_dist
#' ComplexHeatmap::Heatmap(mat_prob_dist$zscore_matrix,cluster_columns = F,cluster_rows = F,show_column_names = F,show_row_names = F,column_title = "z scores") + 
#' ComplexHeatmap::Heatmap(mat_prob_dist$percentile_matrix,cluster_columns = F,cluster_rows = F,show_column_names = F,show_row_names = F,column_title = "percentile scores") + 
#' ComplexHeatmap::Heatmap(mat_prob_dist$original_matrix,cluster_columns = F,cluster_rows = F,show_column_names = F,show_row_names = F,column_title = "original data")
calcCNVKernelProbDist<-function(submatrix=NULL,win=5)
{
  #win<-5
  k = list(matrix=matrix(1/(win^2),nrow=win,ncol=win),kernel='custom')
  class(k)='convKern'
  submatrix_win_avg = spatialfil::applyFilter(submatrix,k)
  diag_avg_matrix<-matrix(0,ncol=ncol(submatrix),nrow=nrow(submatrix))
  diag_sd_matrix<-matrix(0,ncol=ncol(submatrix),nrow=nrow(submatrix))
  for (y in c(nrow(submatrix),1)) 
  {
    coladjustments<-foreach::foreach(x=1:ncol(submatrix),.export=c("y","submatrix","diag_avg_matrix","diag_sd_matrix"),.combine="rbind" ) %do%
    {
      off_diag_for_point<-submatrix[row(submatrix)==col(submatrix)-(y-x)]
      diag_avg<-mean(off_diag_for_point)
      diag_sd<-sd(off_diag_for_point)
      diag_avg_matrix[row(diag_avg_matrix)==col(diag_avg_matrix)-(y-x)]<-diag_avg
      diag_sd_matrix[row(diag_sd_matrix)==col(diag_sd_matrix)-(y-x)]<-diag_sd
      #if(x==319){Heatmap(diag_sd_matrix,cluster_columns = F,cluster_rows=F,show_row_names = F,show_column_names = F)}
      diff<-(y-x)
      debuginfo<-c(diff,y,x)
      names(debuginfo)<-c("diff","y","x")
      debuginfo
    }
  }
  zscore_matrix<-(submatrix_win_avg-diag_avg_matrix)/(diag_sd_matrix)
  zscore_matrix[1,ncol(zscore_matrix)]<-0
  zscore_matrix[ncol(zscore_matrix),1]<-0
  percentile_matrix<-pnorm((submatrix_win_avg-diag_avg_matrix)/diag_sd_matrix)
  output_list<-list(zscore_matrix,percentile_matrix,submatrix)
  names(output_list)<-c("zscore_matrix","percentile_matrix","original_matrix")
  return(output_list)
}
