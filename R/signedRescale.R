#' Rescale positive and negative data, preserving sign information.
#'
#' Performs a signed rescale on the data, shrinking the negative and positive ranges into the [0,1] space, such that negative is always less than 0.5 and positive is always greater.
#' @keywords signed rescale positive negative matrix
#' @import
#' @param matrix A matrix to be transformed
#' @return transformedmatrix A transformed matrix.
#' @export
signedRescale<-function(matrix)
{
  #matrix<-as.matrix(matrix)
  transformedmatrix<-as.matrix(matrix)
  transformedmatrix[transformedmatrix==0]<-0.5+.Machine$double.eps
  transformedmatrix[transformedmatrix>0]<-((transformedmatrix[transformedmatrix>0]/(max(transformedmatrix[transformedmatrix>0])*2))+0.5)
  transformedmatrix[matrix<0]<-((transformedmatrix[transformedmatrix<0]/(min(transformedmatrix[transformedmatrix<0])*2)))
  transformedmatrix[transformedmatrix<0.5 & transformedmatrix>0]<-abs(0.5-transformedmatrix[transformedmatrix<0.5 & transformedmatrix>0])
  transformedmatrix[transformedmatrix==0]<-0.5-.Machine$double.eps
  return(transformedmatrix)
}