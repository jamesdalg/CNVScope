
#' Get the genes in the genomic ranges indicated by the row and column labels.
#'
#' Gets the genes in the ranges within each cell of the matrix.
#' @keywords genomic matrix
#' @import GenomicRanges biomaRt foreach doMC
#' @param genomic_matrix A matrix with row and column names of the format chr1_100_200 (chr,start,end)
#' @param prot_only Inlcude only the protein coding genes from ensembl?
#' @param sequential Turn off parallelism with DoMC?
#' @return concatenated_gene_matrix A matrix with row and column genes
#' @export
getAnnotationMatrix<-function(genomic_matrix,prot_only=T,sequential=F)
{
  if(!exists("grch37")){
  grch37 = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  }
  if(!exists("ensembl_gene_tx_table") | !exists("ensembl_gen_tx_table$gene_biotype"))
{  ensembl_gene_tx_table <- biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id","chromosome_name","transcript_start","transcript_end","start_position","end_position", "strand", "percentage_gene_gc_content","external_gene_name","gene_biotype"),
                                 mart = grch37)
  }
  ensembl_gene_tx_table_prot<-ensembl_gene_tx_table[ensembl_gene_tx_table$gene_biotype=="protein_coding",]
  ensembl_gene_gr<-GenomicRanges::GRanges(seqnames = paste0("chr",ensembl_gene_tx_table$chromosome_name),ranges = IRanges(start = ensembl_gene_tx_table$start_position,end=ensembl_gene_tx_table$end_position),strand = ensembl_gene_tx_table$strand,...=ensembl_gene_tx_table)
  
  ensembl_gene_gr_prot<-GenomicRanges::GRanges(seqnames = paste0("chr",ensembl_gene_tx_table_prot$chromosome_name),
                                               ranges = IRanges(start = ensembl_gene_tx_table_prot$start_position,
                                                                end=ensembl_gene_tx_table_prot$end_position),
                                               strand = ensembl_gene_tx_table_prot$strand,
                                               ...=ensembl_gene_tx_table_prot)
  rownames_gr_genomic_matrix<-underscored_pos_to_GRanges(rownames(genomic_matrix))
  colnames_gr_genomic_matrix<-underscored_pos_to_GRanges(colnames(genomic_matrix))
  if(sequential){registerDoSEQ()} else {registerDoMC()}
  if(prot_only)
  {
    print("prot_only")
    row_gene_strings_genomic_matrix<-foreach(i=1:length(rownames_gr_genomic_matrix),.inorder=T) %dopar% {
      print(i)
      outputstring<-paste(
        unique(gsub("\\..*[[:space:]]","",unique(mcols(subsetByOverlaps(ensembl_gene_gr_prot,rownames_gr_genomic_matrix[i]))$....external_gene_name)))
        ,collapse=" ")
      if(is.null(outputstring) | anyNA(outputstring) | length(outputstring)==0) {outputstring<-""}
      outputstring
    }
    col_gene_strings_genomic_matrix<-foreach(i=1:length(colnames_gr_genomic_matrix),.inorder=T) %dopar% {
      print(i)
      outputstring<-paste(
        unique(gsub("\\..*[[:space:]]","",unique(mcols(subsetByOverlaps(ensembl_gene_gr_prot,colnames_gr_genomic_matrix[i]))$....external_gene_name)))
        ,collapse=" ")
      if(is.null(outputstring) | anyNA(outputstring) | length(outputstring)==0) {outputstring<-""}
      outputstring
    }   
  } else {
    print("all_genes")
    row_gene_strings_genomic_matrix<-foreach(i=1:length(rownames_gr_genomic_matrix),.inorder=T) %dopar% {
      print(i)
      outputstring<-paste(
        unique(gsub("\\..*[[:space:]]","",unique(mcols(subsetByOverlaps(ensembl_gene_gr,rownames_gr_genomic_matrix[i]))$....external_gene_name)))
        ,collapse=" ")
      if(is.null(outputstring) | anyNA(outputstring) | length(outputstring)==0) {outputstring<-""}
      outputstring
    }
    col_gene_strings_genomic_matrix<-foreach(i=1:length(colnames_gr_genomic_matrix),.inorder=T) %dopar% {
      print(i)
      outputstring<-paste(
        unique(gsub("\\..*[[:space:]]","",unique(mcols(subsetByOverlaps(ensembl_gene_gr,colnames_gr_genomic_matrix[i]))$....external_gene_name)))
        ,collapse=" ")
      if(is.null(outputstring) | anyNA(outputstring) | length(outputstring)==0) {outputstring<-""}
      outputstring
    }
  }
  col_gene_strings_matrix_genomic_matrix_alt<-matrix(rep((unlist(col_gene_strings_genomic_matrix)),nrow(genomic_matrix)),ncol=ncol(genomic_matrix),nrow=nrow(genomic_matrix),byrow = T)
  row_gene_strings_matrix_genomic_matrix<-matrix(rep(unlist(row_gene_strings_genomic_matrix),ncol(genomic_matrix)),ncol=ncol(genomic_matrix),nrow=nrow(genomic_matrix)) #essential
  concatenated_gene_matrix<-matrix(
    paste0("row_genes:",row_gene_strings_matrix_genomic_matrix ,"\ncol genes:",col_gene_strings_matrix_genomic_matrix_alt,"\noriginal value:",as.matrix(genomic_matrix))
    ,ncol=ncol(col_gene_strings_matrix_genomic_matrix_alt),
    nrow=nrow(row_gene_strings_matrix_genomic_matrix)) 
  return(concatenated_gene_matrix)
  
}
#genomic_matrix=genomic_matrix
#prot_only=T
#sequential=F
# subset(expressiondata,strsplit(col_gene_strings_matrix_genomic_matrix_alt[1],split=" ")[[1]]
# expression_data<-data.table::fread(paste0(groupdir,"dalgleishjl/hicnv/cor_OS_discovery_exp-cgh.txt"))
#i<-1
#col_gene_strings_matrix_genomic_matrix_alt_split<-sapply(col_gene_strings_matrix_genomic_matrix_alt,function(x) strsplit(x,split=" "))
#single_point_expression<-expression_data[expression_data$gene %in% unlist(col_gene_strings_matrix_genomic_matrix_alt_split[i]) & !(expression_data$chrom %in% c("chrY,chrM")),]
