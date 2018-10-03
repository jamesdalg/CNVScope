#'  Create chromosomal interaction matrices for HiCNV shiny application.
#'
#' Takes a linear regression matrix and sets infinites to a finite value, and changes the sign to match the sign of the correlation for each value.
#' @keywords Interaction matrix
#' @import GenomicInteractions reshape2 magrittr foreach doParallel
#' @importFrom biomaRt useMart getBM
#' @param whole_genome_mat The matrix containing all of the data, from which the individual matrices will be split.
#' @param output_dir the folder where the matrices in RData format, will be written.
#' @param prefix filename prefix for individual matrices. Default: "nbl_"
#' @examples 
#' #examples for this function would be too large to 
#' #include and should be run on an HPC machine node.
#' #illustration of this process is shown clearly in 
#' #the vignette and can be done if a user properly
#' #follows the instructions.
#' # The function is intended to be run on a whole interactome matrix (chr1-X).
#' @return The list of files already written to disk, with full filenames and paths.
#' @export
createChromosomalMatrixSet<-function(whole_genome_mat,output_dir,prefix="nbl_")
{
  original_dir<-getwd()
  if(!dir.exists(output_dir)){dir.create(output_dir)}
  setwd(output_dir)
  chromosomes<-paste0("chr",c(seq(1:22),"X"),"_")
  chrom.pairs<-expand.grid(1:length(chromosomes),1:length(chromosomes))
  grch37 = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  ensembl_gene_tx_table <- biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id","chromosome_name","transcript_start","transcript_end","start_position","end_position", "strand", "percentage_gene_gc_content","external_gene_name","gene_biotype"),
                                          mart = grch37)
  
  foreach(q=rev(1:nrow(chrom.pairs)),.export=ls(),.errorhandling = "pass") %dopar%
  {
    if(!file.exists(paste0(chromosomes[chrom.pairs[q,1]],chromosomes[chrom.pairs[q,2]],prefix,"sample_matched","_unrescaled",".RData")))
    { print(paste0("line of incomplete matrix:",q,chromosomes[chrom.pairs[q,1]],chromosomes[chrom.pairs[q,2]])) 
      full_submatrix_dimensions<-c(length(grep(chromosomes[chrom.pairs[q,1]],rownames(whole_genome_mat))),length(grep(chromosomes[chrom.pairs[q,2]],colnames(whole_genome_mat))))
      ggplotmatrix<-writeAsymmetricMeltedChromosomalMatrixToDisk(whole_genome_matrix = whole_genome_mat,chrom1 = chrom.pairs[q,1],chrom2 = chrom.pairs[q,2],extra_data_matrix = NULL,transpose=T,sequential=T,debug=F,desired_range_start = min(full_submatrix_dimensions),desired_range_end = max(full_submatrix_dimensions),rescale = F,saveToDisk = F) #,flip_row_col_genes = T
      ggplotmatrix<-ggplotmatrix[,c("Var1","Var2","value","value1")]
      ggplotmatrix$Var1<-as.character(ggplotmatrix$Var1)
      ggplotmatrix$Var2<-as.character(ggplotmatrix$Var2)
      ggplotmatrix$value<-as.numeric(unlist(ggplotmatrix$value))
      ggplotmatrix$value1<-as.character(ggplotmatrix$value1)
      ggplotmatrix$orig_value<-reshape2::colsplit(ggplotmatrix$value1,"value:",c("junk","orig_value"))$orig_value
      ggplotmatrix$value<-as.numeric(ggplotmatrix$orig_value)
      save("ggplotmatrix",file=paste0(chromosomes[chrom.pairs[q,1]],chromosomes[chrom.pairs[q,2]],prefix,"sample_matched","_unrescaled",".RData"))
    }
  } 
  outputfilelist<-list.files(pattern=utils::glob2rx("*chr*.RData"))
  setwd(original_dir)
  return(outputfilelist)
}
