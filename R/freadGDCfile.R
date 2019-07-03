#' Read GDC segmentation datafile for low-pass sequencing data.
#'
#' Reads a GDC segmetnation file and extract the segmetnation data.
#' @keywords read file
#' @import magrittr
#' @importFrom data.table fread
#' @importFrom reshape2 colsplit
#' @importFrom dplyr bind_cols
#' @importFrom tibble as.tibble
#' @param file GDC file to be read
#' @param fread_skip The number of metadata lines to be skipped(typically 14)
#' @param format The format of the files (TCGA or TARGET)
#' @references https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
#' @return input_tsv_with_sample_info A data frame containing the sample information extracted
#'  from the filename, including sample name & comparison type.
#' @examples 
#' freadGDCfile(file =
#' system.file("extdata","somaticCnvSegmentsDiploidBeta_TARGET-30-PANRVJ_NormalVsPrimary.tsv",
#' package = "CNVScope"))
#' @export
globalVariables(c('....uuid','barcode1','barcode2','current_gr.....Segment_Mean','fn', 'sep', 'uuid'))
#global variables calls were put in to make it pass CRAN checks. Feel free to disable as needed.  
freadGDCfile<-function(file,fread_skip=NULL, format = "TARGET") {

if(format=="TARGET"){
if(is.null(fread_skip)){  fread_skip=14}
input_tsv<-data.table::fread(file,skip=fread_skip)
sample_info_colsplit<-reshape2::colsplit(basename(file),"_|-|\\.",c("pre","project","num","sample","comparison","fn_ext"))
input_tsv_with_sample_info<-dplyr::bind_cols(input_tsv,sample_info_colsplit[rep(1,nrow(input_tsv)),])
if(length(na.omit(unlist(sample_info_colsplit)))!=6){return(NULL)}
}
if(format=="TCGA")
{
  if(is.null(fread_skip)){  fread_skip=0}
  input_tsv_with_sample_info<-data.table::fread(file,skip=fread_skip) %>% 
  dplyr::mutate(fn = basename(file),sep="--") %>% 
 tidyr::separate(fn,sep="---", into = c("barcode1","barcode2","dataformat")) %>% 
 tidyr::separate(barcode1,sep="-|_", into = c("project1","tss1","participant1","sample1","portion_analyte1","plate1","center1","barcode_id1"),remove=FALSE) %>%
 tidyr::separate(barcode2,sep="-|_", into = c("project2","tss2","participant2","sample2","portion_analyte2","plate2","center2","barcode_id2"),remove=FALSE) %>%
    dplyr::mutate( uuid = dirname(file) %>% tibble::as.tibble() %>% tidyr::separate(value, sep="/",into=c("dir","uuid")) %>% dplyr::pull(uuid)) %>% 
    dplyr::select( -sep)
}
return(input_tsv_with_sample_info)
}
#for reference, see TCGA barcode documentation.