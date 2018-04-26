#' Assign GenomicInteractions to a predefined series of bins for row and column, corresponding to a genomic matrix.
#'
#' This function allows the user to assign a set of genomicinteractions to a pre-existing matrix with known dimensions and column/row names. It finds the row/column index of each point and produces a merged dataframe with the original annotation columns that correspond to each bin in the matrix, with appropriate labels & indexes.
#' @param gint A GenomicInteractions object needing to be binned.
#' @param whole_genome_matrix A matrix with underscored positions for column and rownames e.g. chr1_1_5000,chr1_5001_10000.  If this is provided, it will override rown/column names and GRanges objects.
#' @param rownames_gr A Genomic Ranges object created from the whole genome matrix row names in chr_start_end format, e.g. chr1_1_5000. No effect if whole_genome_mattrix is specified.
#' @param colnames_gr A Genomic Ranges object created from the whole genome matrix column names in chr_start_end format. No effect if whole_genome_mattrix is specified.
#' @param rownames_mat The row names of the whole_genome_matrix in chr_start_end format.
#' @param colnames_mat The column names of the whole_genome_matrix in chr_start_end format.
#' @keywords GenomicInteractions bin matrix colnames rownames binning bin
#' @import foreach doMC GenomicFeatures data.table
#' @export
#' @examples
rebinGenomicInteractions<-function(gint=NULL,whole_genome_matrix=NULL,rownames_gr=NULL,colnames_gr=NULL,rownames_mat=NULL,colnames_mat=NULL)
{
  if(is.null(gint)){return("No GenomicInteractions to rebin!")}
  if(!is.null(whole_genome_matrix)){
    rownames_mat<-rownames(whole_genome_matrix)
    colnames_mat<-rownames(whole_genome_matrix)
  }
  output<-foreach(i=1:length(gint),.inorder = T,.combine="rbind",.errorhandling = "pass",.export = ls()) %dopar% #length(breakpoint_gint_full)length(breakpoint_gint_full)
  {
    #current_int_df<-as.data.table(gint[i]) #
    current_int_df<-as.data.frame(cbind(as.data.frame(anchorOne(gint)[i]),as.data.frame(anchorTwo(gint)[i]),as.data.frame(mcols(gint[i]))))
    print(paste0(i/length(gint)*100,"%"))
    row_bin_index<-findOverlaps(rownames_gr,anchorOne(gint[i]))@from
    col_bin_index<-findOverlaps(colnames_gr,anchorTwo(gint[i]))@from
    if(length(row_bin_index)==0 | length(col_bin_index)==0) {return(NULL)}
    col_bin_label<-colnames_mat[col_bin_index]
    row_bin_label<-rownames_mat[row_bin_index]
    outputline<-c(row_bin_index,col_bin_index,row_bin_label,col_bin_label,sapply(current_int_df,as.character))
    outputnames<-c("row_bin_index","col_bin_index","row_bin_label","col_bin_label",colnames(current_int_df))
    if(length(row_bin_index)==0 | length(col_bin_index)==0) {outputline<-rep("",length(outputnames))}
    names(outputline)<-outputnames
    outputline
  }
  output_df<-as.data.frame(matrix(as.character(output),nrow=nrow(output)))
  colnames(output_df)<-c("row_bin_index","col_bin_index","row_bin_label","col_bin_label",as.character(colnames(as.data.frame(gint[1]))))
  return(output_df)
}