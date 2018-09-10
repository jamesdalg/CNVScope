#' Convert GRanges object to underscord positions.
#'
#' This function converts row or column names (or any character vector of the format) into a GenomicRanges object.
#' @param input_gr A GenomicRanges object
#' @param minusOneToEnd Minus one position to end of each Genomic Range?
#' @keywords Genomic Ranges position
#' @import GenomicFeatures
#' @export
#' @examples
GRanges_to_underscored_pos<-function(input_gr,minusOneToEnd=T)
{
  if(minusOneToEnd){adjustment<-1} else {adjustment=0}
  output_char<-paste0(seqnames(input_gr),"_",input_gr@ranges@start,"_",input_gr@ranges@start+input_gr@ranges@width-adjustment)
  return(output_char)
}