#' Convert GRanges object to underscord positions.
#'
#' This function converts row or column names (or any character vector of the format) into a GenomicRanges object.
#' @param input_gr A GenomicRanges object
#' @param minusOneToEnd Minus one position to end of each Genomic Range?
#' @keywords Genomic Ranges position
#' @importFrom GenomicRanges seqnames GRanges
#' @examples 
#' load(system.file("extdata","nbl_result_matrix_sign_small.rda",package = "HiCNV")) 
#' col_gr<-underscored_pos_to_GRanges(colnames(nbl_result_matrix_sign_small))
#' GRanges_to_underscored_pos(col_gr)
#' @export
GRanges_to_underscored_pos<-function(input_gr,minusOneToEnd=T)
{
  if(minusOneToEnd){adjustment<-1} else {adjustment=0}
  output_char<-paste0(seqnames(input_gr),"_",input_gr@ranges@start,"_",input_gr@ranges@start+input_gr@ranges@width-adjustment)
  return(output_char)
}