#' Import a breakpoint BED file.
#'
#' Imports a BED file with breakpoints or other interactions, in a dual position format.
#'
#' @name importBreakpointBed
#' @keywords bed
#' @import GenomicInteractions
#' @importFrom rtracklayer import.bed
#' @importFrom reshape2 colsplit
#' @param breakpoint_fn the filename of the breakpoint bed file
#' @return a Genomic Interactions Object
#' @examples 
#' importBreakpointBed(breakpoint_fn = system.file("extdata",
#' "sample_breakpoints.bed",package = "CNVScope"))
#' @export
globalVariables(c("mcols","mcols<-"))
  importBreakpointBed<-function(breakpoint_fn)
  {
    imported_bed<-rtracklayer::import.bed(breakpoint_fn)
    colsplit_locations<-reshape2::colsplit(string = GenomicRanges::mcols(imported_bed)$name,pattern = "_R_",names=c("otherdata","pos"))[,2]
    colsplit_locations_gr<-GRanges(colsplit_locations)
    colsplit_id<-reshape2::colsplit(string=reshape2::colsplit(string = GenomicRanges::mcols(imported_bed)$name,pattern = "Id:",names=c("otherdata","id_containing"))[,2],pattern="_",names=c("id","other"))[,1]
    gint<-GenomicInteractions(imported_bed,colsplit_locations_gr,... = colsplit_id)
    colsplit_bed<-reshape2::colsplit(string = GenomicRanges::mcols(imported_bed)$name,pattern = "_",names=c("otherdata","id_containing"))[,1]
    S4Vectors::mcols(gint)$id<-colsplit_id
    S4Vectors::mcols(gint)$bed<-colsplit_bed 
    return(gint)
  }
