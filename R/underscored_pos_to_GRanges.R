#' Convert coordinates in underscored format to a GRanges object.
#'
#' This function creates a new GRanges object from a character vector of coordinates in the form "chr1_0_5000" and creates a GRanges object from them.
#' @keywords CNV GRanges Genomic Ranges position
#' @param underscored_positions A vector of positions of the form c("chr1_0_5000","chr1_7500_10000","chr1_10000_15000")
#' @param extended_data Optional metadata columns. These columns cannot be named "start", "end", "width", or "element". Passed to GRanges object as ...
#' @param zeroToOneBasedStart Converts a set of underscored positions that begin with zero to GRanges where the lowest positional value on a chromosome is 1. Essentially adds 1 to start
#' @param zeroToOneBasedEnd Adds 1 to the end of the underscored positions
#' @return A GRanges object
#' @examples
#' load(system.file("extdata","nbl_result_matrix_sign_small.rda",package = "CNVScope"))
#' underscored_pos_to_GRanges(colnames(nbl_result_matrix_sign_small))
#' @export
underscored_pos_to_GRanges<-function(underscored_positions=NULL,extended_data=NULL,zeroToOneBasedStart=T,zeroToOneBasedEnd=F)
{
  #importFrom GenomicRanges GRanges
  #importFrom IRanges IRanges
  #importFrom plyr ldply
  GenomicRanges::GRanges(seqnames = as.character(plyr::ldply(sapply(underscored_positions,
                                                      function(x) strsplit(as.character(x),"_",fixed=T)
),rbind)["1"][,1]),IRanges::IRanges(start = 
                             as.numeric(as.character(plyr::ldply(sapply(underscored_positions,
                                                                        function(x) strsplit(as.character(x),"_",fixed=T)
                             ),rbind)["2"][,1])) +   as.numeric(zeroToOneBasedStart) ,end = as.numeric(as.character(plyr::ldply(sapply(underscored_positions,
                                                                                                  function(x) strsplit(as.character(x),"_",fixed=T)
                             ),rbind)["3"][,1])) +   as.numeric(zeroToOneBasedEnd)),... = extended_data)
}
