#' Average edges of a matrix to facilitate downsampling.
#'
#' Averages the columns and rows of a matrix by a certain amount.
#' @keywords rescale downsample average edges matrix
#' @importFrom reshape2 colsplit
#' @importFrom Matrix colMeans rowMeans
#' @param unchangedmatrix A matrix to have edges averaged with genomic coordinates in the form chr1_50_100 set as the column and row names.
#' @param nedges The number of edges to be averaged
#' @param dimension Selectively averages edges in one dimension. Performs symmetric edge averaging by default.
#' @return averaged_matrix A matrix with edges averaged, which may be more amenable to downsampling
#' @export
averageMatrixEdges<-function(unchangedmatrix,nedges=1,dimension=c("row","column"))
{
  #dimension<-gsub("col","column",dimension)
  if(!(length(intersect(dimension,"row"))==1 | length(intersect(dimension,"column"))==1 )) { 
    errormsg<-paste0("Invalid dimension specification:\'",dimension,"\' Valid options are \'column\' and \'row\'")
    #if(length)
    print(paste0(dimension))
    print(errormsg)
    stop()
    return(errormsg)}
  if("row" %in% dimension)
  {
    #dim(unchangedmatrix[(nrow(unchangedmatrix)-nedges):nrow(unchangedmatrix),])
    #length(colMeans(unchangedmatrix[(nrow(unchangedmatrix)-nedges):nrow(unchangedmatrix),]))
    averaged_row<-(Matrix::colMeans(unchangedmatrix[(nrow(unchangedmatrix)-nedges):nrow(unchangedmatrix),]))
    averaged_rownames_df<-reshape2::colsplit(string = rownames(unchangedmatrix)[(nrow(unchangedmatrix)-nedges):nrow(unchangedmatrix)],pattern = "_",names = c("chrom","start","end"))
    #<-(nrow(unchangedmatrix)-nedges):nrow(unchangedmatrix)
    #averaged_matrix<-unchangedmatrix[1:nrow(unchangedmatrix)-nedges-1),]
    #dim(unchangedmatrix[1:(nrow(unchangedmatrix)-nedges-1),])
    averaged_matrix<-rbind(unchangedmatrix[1:(nrow(unchangedmatrix)-nedges-1),],averaged_row)
    rownames(averaged_matrix)[nrow(averaged_matrix)]<-paste(c(as.character(averaged_rownames_df[1,c("chrom","start")]),as.character(averaged_rownames_df[nrow(averaged_rownames_df),c("end")])),collapse = "_")
    #dim(averaged_matrix)
    original_matrix<-unchangedmatrix
    unchangedmatrix<-averaged_matrix
    if(!("column" %in% dimension))
    {
      return(averaged_matrix)
    }
  }
  
  if("column" %in% dimension)
  {
    #dim(unchangedmatrix)
    averaged_column<-(Matrix::rowMeans(unchangedmatrix[,(ncol(unchangedmatrix)-nedges):ncol(unchangedmatrix)]))
    averaged_colnames_df<-reshape2::colsplit(string = colnames(unchangedmatrix)[(ncol(unchangedmatrix)-nedges):ncol(unchangedmatrix)],pattern = "_",names = c("chrom","start","end"))
    averaged_matrix<-cbind(unchangedmatrix[,1:(ncol(unchangedmatrix)-nedges-1)],averaged_column)
    colnames(averaged_matrix)[ncol(averaged_matrix)]<-paste(c(as.character(averaged_colnames_df[1,c("chrom","start")]),as.character(averaged_colnames_df[nrow(averaged_colnames_df),c("end")])),collapse = "_")
    return(averaged_matrix)
  }
}