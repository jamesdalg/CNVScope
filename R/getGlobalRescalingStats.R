#'  Calculate several base statistics for color rescaling.
#'
#' calculates several statistics from a large matrix that can then be applied to smaller submatrices without needing to load the entire matrix into memmory
#' @keywords rescale color stats
#' @import stats
#' @param whole_matrix the whole matrix to get stats for.
#' @param saveToDisk Save the statistics to disk as an RDS file in the local directory?
#' @param output_fn the name of the output file.
#' @return A list of the output statistics, including:
#' the global min, max, length, sigma (matrix variance), pos_sigma (variance of the positive values), neg_sigma(variance of the negative values), global mean (global_mu),
#'  est_max_cap (global_mu+global_sigma_pos*2), as well as the number of rows and columns of the matrix.
#' @export
getGlobalRescalingStats<-function(whole_matrix,saveToDisk=F,output_fn=NULL)
{
  if(is.null(output_fn)){output_fn<-"whole_matrix_stats.rds"}
  whole_matrix<-as.matrix(whole_matrix)
global_max<-max(whole_matrix)
global_min<-min(whole_matrix)
global_length<-length(whole_matrix)
global_sigma<-sqrt(var(whole_matrix)*(global_length)/(global_length-1))
global_sigma_pos<-sqrt(var(whole_matrix[whole_matrix>0])*(global_length)/(global_length-1))
global_sigma_neg<-sqrt(var(whole_matrix[whole_matrix<0])*(global_length)/(global_length-1))
global_mu<-mean(whole_matrix)
est_max_cap<-global_mu+global_sigma_pos*2
nrow_mat<-nrow(whole_matrix)
ncol_mat<-ncol(whole_matrix)
output<-list(global_max,global_min,global_sigma_pos,global_sigma_neg,global_mu,est_max_cap,nrow_mat,ncol_mat,global_length)
if(saveToDisk){saveRDS(object = "output",file=paste0(output_fn))}

return(output)
}
#global_stats<-getGlobalRescalingStats(all_conc_cleaned_common_coords_linreg)