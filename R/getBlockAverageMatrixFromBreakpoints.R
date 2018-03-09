#' Calculate block averages and areas in a matrix given breakpoints.
#'
#' This function produces several matrix outputs of averages and areas of matrix blocks, given a pair of vectors for breakpoints.
#' @keywords CNV kernel probability distribution concordance fast
#' @import ComplexHeatmap foreach doMC spatialfil 
#' @param whole_matrix the large, whole matrix from which blocks are taken
#' @param breakpoints_col An integer list of column breakpoints, including 1 and the number of columns in the whole matrix.
#' @param breakpoints_row An integer list of row breakpoints, including 1 and the number of rows in the whole matrix.
#' @return An output list of the following:
#' @return blockaverages_reformatted_by_index  a matrix of the block averages and areas, in long format, with indexes used to generate the averages.
#' @return blockaverages_reformatted_by_label a matrix of the block averages and areas, in long format, with labels of the indexes used to generate the averages.
#' @return blockaverages_matrix_idx_area a matrix of the block areas, with indexes based on the original row/col index used to generate the data.
#' @return blockaverages_matrix_idx_avg a matrix of the block averages, with indexes based on the original row/col index used to generate the data.
#' @return blockaverages_matrix_label_area a matrix of the block areas, with indexes based on the original row/col label used to generate the data. 
#' @return blockaverages_matrix_label_avg a matrix of the block averages, with indexes based on the original row/col label used to generate the data. 
#' @examples
#' set.seed(303)
#' mat<-matrix(data=runif(n = 25),nrow=5,ncol=5,dimnames = list(c("chr1_0_5000","chr1_5000_10000","chr1_10000_15000","chr1_15000_20000","chr1_20000_25000"),c("chr1_0_5000","chr1_5000_10000","chr1_10000_15000","chr1_15000_20000","chr1_20000_25000")))
#' breakpoints_col<-c(1,2,4,5)
#' breakpoints_row<-c(1,2,4,5)
#' HiCNV::getBlockAverageMatrixFromBreakpoints(whole_matrix=mat,breakpoints_col=breakpoints_col,breakpoints_row=breakpoints_row)
#' @export
getBlockAverageMatrixFromBreakpoints<-function(whole_matrix,breakpoints_col,breakpoints_row,outputs=c("blockaverages_reformatted_by_index","blockaverages_reformatted_by_label","blockaverages_matrix_idx_area","blockaverages_matrix_idx_avg","blockaverages_matrix_label_avg","blockaverages_matrix_label_area"))
{
  breakpoints_col<-as.integer(unique(
    gsub("^0$",1,c(0,breakpoints_col,ncol(whole_matrix)))
  ))
  breakpoints_row<-as.integer(unique(
    gsub("^0$",1,c(0,breakpoints_row,nrow(whole_matrix)))
  ))
  
  blockaverages<-foreach(j=1:(length(breakpoints_col)-1),.combine="rbind",.inorder=T) %do%
  {
    
    foreach(i=1:(length(breakpoints_row)-1),.combine="rbind",.inorder=T) %dopar%
    {
      #whole_matrix[i,j]
      print(paste0("i",i,"j",j))
      t(as.data.frame(c(
        breakpoints_row[i],breakpoints_col[j],mean(as.numeric(unlist(whole_matrix[breakpoints_row[i]:breakpoints_row[i+1],breakpoints_col[j]:breakpoints_col[j+1]]))),
        as.numeric(abs(breakpoints_row[i]-breakpoints_row[i+1])*abs(breakpoints_col[j]-breakpoints_col[j+1]))
      )))
    }
    
  }
  #blockaverages_reformatted<-matrix(as.numeric(blockaverages),ncol=3,nrow=((length(rowindices_subset)-1)*3))
  colnames(blockaverages)<-c("rowindex","colindex","average","area")
  rownames(blockaverages)<-NULL
  #blockaverages[blockaverages[,"area"]==0,]<-NULL
  blockaverages_reformatted_by_index<-blockaverages[(blockaverages[,"area"]!=0),]
  
  blockaverages_reformatted_by_label<-foreach(i=1:nrow(blockaverages_reformatted_by_index),.combine="rbind") %do%
  {
    if(nrow(blockaverages_reformatted_by_index)!=i)
    {outputline<-c(
      colnames(whole_matrix)[as.integer(gsub(0,1,blockaverages[i,1]))],rownames(whole_matrix)[as.integer(gsub(0,1,blockaverages[i,2]))],
      colnames(whole_matrix)[as.integer(gsub(0,1,blockaverages[i+1,1]))],rownames(whole_matrix)[as.integer(gsub(0,1,(blockaverages[i+1,2])))],blockaverages_reformatted_by_index[i,3],blockaverages_reformatted_by_index[i,"area"])}
  }
  colnames(blockaverages_reformatted_by_label)<-c("rowstart","rowend","colstart","colend","average","area")
  #still need to convert this into a matrix
  blockaverages_reformatted_by_index_df<-as.data.frame(blockaverages_reformatted_by_index)
  blockaverages_reformatted_by_label_df<-as.data.frame(blockaverages_reformatted_by_label)
  blockaverages_matrix_idx_avg<-dcast(blockaverages_reformatted_by_index_df,rowindex ~ colindex,value.var = "average")
  blockaverages_matrix_idx_area<-dcast(blockaverages_reformatted_by_index_df,rowindex ~ colindex,value.var = "area")
  blockaverages_matrix_label_avg<-blockaverages_matrix_idx_avg
  blockaverages_matrix_label_area<-blockaverages_matrix_idx_area
  #breakpoints_col_zero_corr<-gsub(0,1,breakpoints_col)
  #breakpoints_row_zero_corr<-gsub(0,1,breakpoints_row)
  colnames(blockaverages_matrix_label_avg)[2:length(colnames(blockaverages_matrix_label_avg))]<-colnames(whole_matrix)[as.integer(colnames(blockaverages_matrix_label_avg)[2:length(colnames(blockaverages_matrix_label_avg))])]
  rownames(blockaverages_matrix_label_avg)<-rownames(whole_matrix)[as.integer(rownames(blockaverages_matrix_label_avg))]
  colnames(blockaverages_matrix_label_area)[2:length(colnames(blockaverages_matrix_label_area))]<-colnames(whole_matrix)[as.integer(colnames(blockaverages_matrix_label_area)[2:length(colnames(blockaverages_matrix_label_area))])]
  rownames(blockaverages_matrix_label_area)<-rownames(whole_matrix)[as.integer(rownames(blockaverages_matrix_label_area))]
  #  blockaverages_matrix_label<-dcast(blockaverages_reformatted_by_label_df[,c("rowstart","colstart","average")],rowstart ~ colstart,value.var = "average",fun.aggregate = mean)
  outputlist<-list(blockaverages_reformatted_by_index,blockaverages_reformatted_by_label,blockaverages_matrix_idx_area,blockaverages_matrix_idx_avg,blockaverages_matrix_label_avg,blockaverages_matrix_label_area)
  names(outputlist)<-c("blockaverages_reformatted_by_index","blockaverages_reformatted_by_label","blockaverages_matrix_idx_area","blockaverages_matrix_idx_avg","blockaverages_matrix_label_avg","blockaverages_matrix_label_area")
  return(outputlist)
}
