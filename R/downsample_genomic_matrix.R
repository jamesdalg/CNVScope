#' Rescale positive and negative data, preserving sign information.
#'
#' Downsamples a matrix by a specified factor.
#' @name downsample_genomic_matrix
#' @keywords signed rescale positive negative matrix
#' @import reshape2 OpenImageR foreach
#' @param whole_matrix A matrix to be downsampled, on a single chromosome
#' @param downsamplefactor A factor by which to reduce the matrix. Must be something that both the row and columns can be divisible by.
#' @param singlechromosome Single chromosome mode; Multi-chromosome not yet implemented (leave T)
#' @return whole_matrix_dsamp A downsampled matrix.
#' @examples 
#' load(system.file("extdata","nbl_result_matrix_sign_small.rda",package = "CNVScope"))
#' downsample_genomic_matrix(whole_matrix=nbl_result_matrix_sign_small,
#' downsamplefactor=5,singlechromosome=TRUE)
#' @export
#lines to disable spurrious CRAN notes below. i is obviously not global.
globalVariables("i")
dontCheck('i')
downsample_genomic_matrix<-function(whole_matrix,downsamplefactor,singlechromosome=T)
{
  i <- if(exists("i")){get("i")} else {NULL}
  if(singlechromosome)
{  if(nrow(whole_matrix)%%downsamplefactor==0 & ncol(whole_matrix)%%downsamplefactor==0 )
  {
    whole_matrix_dsamp<-down_sample_image(whole_matrix,factor=downsamplefactor,gaussian_blur = T)  
  }  else {(return("downsample not a factor of rows and columns"))}
  #combine labels for the downsampled matrix (for factor 5, take five ranges and combine them into 1.)  
  #essentially take the first part of the matrix (line 1 chr1_1235, then concatenate it to the final point in the range of the last bin).

  downsampled_colnames<-unlist(foreach(i=seq(from=1,to=ncol(whole_matrix),by=downsamplefactor)) %do%
  {
    paste0(paste(reshape2::colsplit(colnames(whole_matrix)[i],"_",names=c("chrom","start","end"))[,c(1,2)],collapse="_"),"_",
           reshape2::colsplit(colnames(whole_matrix)[i+downsamplefactor-1],"_",names=c("chrom","start","end"))[,c(3)])
  })
  downsampled_rownames<-unlist(foreach(i=seq(from=1,to=nrow(whole_matrix),by=downsamplefactor)) %do%
  {
    paste0(paste(reshape2::colsplit(rownames(whole_matrix)[i],"_",names=c("chrom","start","end"))[,c(1,2)],collapse="_"),"_",
           reshape2::colsplit(rownames(whole_matrix)[i+downsamplefactor-1],"_",names=c("chrom","start","end"))[,c(3)])
  })
  colnames(whole_matrix_dsamp)<-downsampled_colnames
  rownames(whole_matrix_dsamp)<-downsampled_rownames
  return(whole_matrix_dsamp)
  }
  if(!singlechromosome)
  {
    return("multi-chromosome not yet implemented")
  }
  }
test_mat<-matrix(runif(81),nrow=9)
colnames(test_mat)<-paste0("chr1_",1:ncol(test_mat))
rownames(test_mat)<-paste0("chr1_",1:nrow(test_mat))
downsample_genomic_matrix(whole_matrix=test_mat,downsamplefactor=3,singlechromosome=T)
