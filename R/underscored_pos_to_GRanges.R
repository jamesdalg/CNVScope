#' Convert coordinates in underscored format to a GRanges object.
#'
#' This function creates a new GRanges object from a character vector of coordinates in the form "chr1_0_5000" and creates a GRanges object from them.
#' @keywords CNV GRanges Genomic Ranges position
#' @import plyr GenomicRanges
#' @param underscored_positions A vector of positions of the form c("chr1_0_5000","chr1_7500_10000","chr1_10000_15000")
#' @param zeroToOneBased Convert a set of underscored positions that begin with zero to GRanges where the lowest positional value on a chromosome is 1.
#' @return A GRanges object
#' @examples
#' HiCNV::underscored_pos_to_GRanges(c("chr1_0_5000","chr1_7500_10000","chr1_10000_15000"))
#' @export
underscored_pos_to_GRanges<-function(underscored_positions=NULL,extended_data=NULL,zeroToOneBased=T)
{
  GRanges(seqnames = as.character(plyr::ldply(sapply(underscored_positions,
                                                      function(x) strsplit(as.character(x),"_",fixed=T)
),rbind)["1"][,1]),IRanges(start = 
                             as.numeric(as.character(plyr::ldply(sapply(underscored_positions,
                                                                        function(x) strsplit(as.character(x),"_",fixed=T)
                             ),rbind)["2"][,1])) +   as.integer(zeroToOneBased),end = as.numeric(as.character(plyr::ldply(sapply(underscored_positions,
                                                                                                  function(x) strsplit(as.character(x),"_",fixed=T)
                             ),rbind)["3"][,1]))),... = extended_data) +   as.integer(zeroToOneBased)
}
