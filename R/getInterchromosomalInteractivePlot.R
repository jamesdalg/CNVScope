#' Create an HTML widget for use in shiny or webshot for a given pair of chromosomes.
#'
#' This function requires a matrix with genomic coordinates in the row and column names, and produces a heatmap with a tooltip
#' @name getInterchromosomalInteractivePlot
#' @keywords CNV heatmap HTML widget data.table readr
#' @importFrom biomaRt getBM useMart
#' @importFrom GenomicRanges GRanges seqnames mcols
#' @importFrom IRanges subsetByOverlaps
#' @importFrom heatmaply heatmaply
#' @importFrom ggplot2 scale_fill_gradient2
#' @param whole_matrix the large, whole genomic matrix from which the submatrix is taken (rows)
#' @param chrom1 The first chromsome used for the map (columns).
#' @param chrom2 The second chromsome used for a map axis.
#' @return An HTML widget.
#' @examples
#' \dontrun{
#' load(system.file("extdata","nbl_result_matrix_sign_small.rda",package = "CNVScope")) 
#' getInterchromosomalInteractivePlot(whole_matrix=nbl_result_matrix_sign_small,chrom1=1,
#' chrom2=1)
#' }
#' @export

globalVariables(c('chromosomes'))
getInterchromosomalInteractivePlot<-function(whole_matrix,chrom1,chrom2)
{
  #if(rownames(whole_matrix)==colnames(whole_matrix))
  chromosomes<-paste0("chr",c(seq(1:22),"X"),"_")
    submatrix<-whole_matrix[grep(chromosomes[chrom1],rownames(whole_matrix)),grep(chromosomes[chrom2],colnames(whole_matrix))]
  
  grch37 = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  
  ensembl_gene_tx_table <- biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id","chromosome_name","transcript_start","transcript_end","start_position","end_position", "strand", "percentage_gene_gc_content","external_gene_name"),
                                          # filters = "ensembl_transcript_id", values = "ENST00000296026",
                                          mart = grch37)
  
  ensembl_gene_gr<-GenomicRanges::GRanges(seqnames = paste0("chr",ensembl_gene_tx_table$chromosome_name),ranges = IRanges(start = ensembl_gene_tx_table$start_position,end=ensembl_gene_tx_table$end_position),strand = ensembl_gene_tx_table$strand,...=ensembl_gene_tx_table)
  
  if(substr(chrom1,start = nchar(chrom1),stop = nchar(chrom1))!="_"){chrom1<-paste0(chrom1,"_")}
  if(substr(chrom2,start = nchar(chrom2),stop = nchar(chrom2))!="_"){chrom1<-paste0(chrom2,"_")}
  
  rownames_gr_submatrix<-underscored_pos_to_GRanges(rownames(submatrix))
  colnames_gr_submatrix<-underscored_pos_to_GRanges(colnames(submatrix))
  row_gene_strings_submatrix<-foreach(i=1:length(rownames_gr_submatrix),.inorder=T) %do% {
    print(i)
    outputstring<-paste(
      unique(gsub("\\..*[[:space:]]","",unique(GenomicRanges::mcols(IRanges::subsetByOverlaps(ensembl_gene_gr,rownames_gr_submatrix[i]))$....external_gene_name)))
      ,collapse=" ")
    if(is.null(outputstring) | anyNA(outputstring) | length(outputstring)==0) {outputstring<-""}
    outputstring
  }
  col_gene_strings_submatrix<-foreach(i=1:length(colnames_gr_submatrix),.inorder=T) %do% {
    print(i)
    outputstring<-paste(
      unique(gsub("\\..*[[:space:]]","",unique(GenomicRanges::mcols(IRanges::subsetByOverlaps(ensembl_gene_gr,colnames_gr_submatrix[i]))$....external_gene_name)))
      ,collapse=" ")
    if(is.null(outputstring) | anyNA(outputstring) | length(outputstring)==0) {outputstring<-""}
    outputstring
  }
  #col_gene_strings_matrix_submatrix<-matrix(rep(unlist(col_gene_strings_submatrix),nrow(submatrix)),ncol=nrow(submatrix),nrow=ncol(submatrix))
  #col_gene_strings_matrix_submatrix_transposed<-t(col_gene_strings_matrix_submatrix)
  col_gene_strings_matrix_submatrix_alt<-matrix(rep((unlist(col_gene_strings_submatrix)),nrow(submatrix)),ncol=ncol(submatrix),nrow=nrow(submatrix),byrow = T) #necessary
  row_gene_strings_matrix_submatrix<-matrix(rep(unlist(row_gene_strings_submatrix),ncol(submatrix)),ncol=ncol(submatrix),nrow=nrow(submatrix)) #necessary
  concatenated_gene_matrix<-matrix(
    paste0("row_genes:",row_gene_strings_matrix_submatrix ,"\ncol genes:",col_gene_strings_matrix_submatrix_alt,"\noriginal value:",as.matrix(submatrix))
    ,ncol=ncol(col_gene_strings_matrix_submatrix_alt),
    nrow=nrow(row_gene_strings_matrix_submatrix)) 
  
  #if(Sys.info()['sysname']=="Windows"){groupdir<-"W:/"} else {groupdir<-"/data/CCRBioinfo/"}
  htmlwidget<-heatmaply(signedRescale(submatrix),Rowv = F,Colv = F,showticklabels=F,custom_hovertext = concatenated_gene_matrix,
                        #file=paste0(groupdir,"dalgleishjl/hicnv/inter/fixed_rescale/largemem2/chr",chromosomes[chrom1],chromosomes[chrom2],"withrow_and_colgenes_genes_fixed_rescale.html"),
                        scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.5, limits = c(0, 1)))
  return(htmlwidget)
}