#' Get Block Indices from an asymmetric (or symmetric) matrix.
#' 
#' This function segments a matrix, including asymmetric matrices using multiple imputation (MI) techniques and a segmentation algorithm to generate breakpoints for column and row.
#' 
#' @keywords HiCseg MI multiple imputation Hi-C CNV breakpoints
#' @import HiCseg blockseg jointseg
#' @param genomicmatrix the large, whole matrix from which blocks are taken
#' @param nb_change_max the maximal number of changepoints, passed to HiCseg (if this algorithm is used). Note: HiCseg doesn't actually obey this limit. Rather, use it as a parameter to increase/decrease segmentation extent.
#' @param distrib Passed to Hicseg_linkC_R, from their documentation: Distribution of the data: "B" is for Negative Binomial distribution, "P" is for the Poisson distribution and "G" is for the Gaussian distribution."
#' @param model Passed on to HiCseg_linkC_R: "Type of model: "D" for block-diagonal and "Dplus" for the extended block-diagonal model."
#' @param MI_strategy strategy to make the matrix temporarily symmetric. "average" adds a number of values equal to the average of the matrix, while copy copies part of the matrix to the shorter side, making a square matrix.
#' @param transpose transpose the matrix and output the breakpoints? Some segmentation algorithms (e.g. HiCseg) produces different results when used against the transposed version of the matrix, as it expects symmetry. This allows the output of additional breakpoints Users can choose to take intersect() or union() on the results to get conserved changepoints or additional changepoints, depending on need.
#' @return An output list of the following:
#' @return breakpoints_col A vector of breakpoints for the columns.
#' @return breakpoints_row A vector of breakpoints for the rows.
#' @return breakpoints_col A vector of breakpoints for the columns on the transposed genomic matrix.
#' @return breakpoints_row A vector of breakpoints for the rows on the transposed genomic matrix.
#' @examples 
#' \dontrun{
#' submatrix_tiny #from chr17 vs chr 6, 5x10.
#' tiny_test<-getAsymmetricBlockIndices(submatrix_tiny)
#' submatrix_wide<-submatrix[1:5,]
#' # > dim(submatrix_wide)
#' # [1]   5 856
#' submatrix_narrow<-submatrix[.1:5]
#' wide_test<-getAsymmetricBlockIndices(submatrix_wide,distrib = "G",model = "Dplus",
#' nb_change_max = 1e4)
#' narrow_test<-getAsymmetricBlockIndices(submatrix_narrow,distrib="P",model="D")
#' random_wide<-matrix(runif(n = 400*200),ncol=400,nrow=200)
#' random_narrow<-matrix(runif(n = 400*200),ncol=200,nrow=400)
#' random_wide_test_avg<-getAsymmetricBlockIndices(random_wide,
#' distrib = "G",model = "Dplus",nb_change_max = 1e4)
#' random_narrow_test_avg<-getAsymmetricBlockIndices(random_narrow,
#' distrib = "G",model = "Dplus",nb_change_max = 1e4)
#' random_wide_test_copy<-getAsymmetricBlockIndices(random_wide,
#' distrib = "G",model = "Dplus",nb_change_max = 1e4,MI_strategy = "copy")
#' random_narrow_test_copy<-getAsymmetricBlockIndices(random_narrow,
#' distrib = "G",model = "Dplus",nb_change_max = 1e4,MI_strategy = "copy")
#' genomicmatrix=random_narrow
#' nb_change_max=100
#' model = "D"
#' distrib = "G"
#' MI_strategy="copy"
#' #question-- does it pick different breakpoints if transposed first? 
#' #Answer: yes, at least in Dplus model.
#' rm(genomicmatrix)
#' rm(model)
#' rm(distrib)
#' rm(MI_strategy)
#' random_wide_test_copy<-getAsymmetricBlockIndices(genomicmatrix = random_wide,
#' distrib = "G",model = "Dplus",nb_change_max = 1e2,MI_strategy = "copy")
#' random_narrow_test_copy<-getAsymmetricBlockIndices(random_narrow,distrib = "G",
#' model = "Dplus",nb_change_max = 1e2,MI_strategy = "copy")
#' random_wide_test_copy_t<-getAsymmetricBlockIndices(genomicmatrix = t(random_wide),
#' distrib = "G",model = "Dplus",nb_change_max = 1e2,MI_strategy = "copy")
#' random_narrow_test_copy_t<-getAsymmetricBlockIndices(genomicmatrix = t(random_narrow),
#' distrib = "G",model = "Dplus",nb_change_max = 1e2,MI_strategy = "copy")
#' length(intersect(random_wide_test_copy$breakpoints_col,
#' random_wide_test_copy_t$breakpoints_row))/length(unique(c(random_wide_test_copy$breakpoints_col,
#' random_wide_test_copy_t$breakpoints_row)))
#' random_wide_test_copy_with_transpose<-getAsymmetricBlockIndices(genomicmatrix = random_wide,
#' distrib = "G",model = "Dplus",nb_change_max = 1e2,MI_strategy = "copy",transpose = T)
#' random_narrow_test_copy_with_transpose<-getAsymmetricBlockIndices(genomicmatrix = random_narrow,
#' distrib = "G",model = "Dplus",nb_change_max = 1e2,MI_strategy = "copy",transpose = T)
#' #passes tests
#' random_narrow_test_copy_with_transpose<-getAsymmetricBlockIndices(genomicmatrix = random_narrow,
#' distrib = "G",model = "Dplus",nb_change_max = 1e2,MI_strategy = "copy",transpose = T)
#' conserved_breakpoints_col<-intersect(random_narrow_test_copy_with_transpose$breakpoints_col,
#' random_narrow_test_copy_with_transpose$t_breakpoints_row)
#' conserved_breakpoints_row<-intersect(random_narrow_test_copy_with_transpose$breakpoints_row,
#' random_narrow_test_copy_with_transpose$t_breakpoints_col)
#' random_wide_test_copy_with_transpose<-getAsymmetricBlockIndices(genomicmatrix = random_wide,
#' distrib = "G",model = "Dplus",nb_change_max = 1e2,MI_strategy = "copy",transpose = T)
#' conserved_breakpoints_col<-intersect(random_wide_test_copy_with_transpose$breakpoints_col,
#' random_wide_test_copy_with_transpose$t_breakpoints_row)
#' conserved_breakpoints_row<-intersect(random_wide_test_copy_with_transpose$breakpoints_row,
#' random_wide_test_copy_with_transpose$t_breakpoints_col)
#' }
#' @export
getAsymmetricBlockIndices<-function(genomicmatrix=NULL,nb_change_max=100,distrib = "G",model = "D",MI_strategy="average",transpose=T)
{
  
  #converting from Asymmetric to symmetric
  if(nrow(genomicmatrix)==ncol(genomicmatrix))
  { 
    hicsegresults<-HiCseg::HiCseg_linkC_R(size_mat=dim(genomicmatrix)[1],nb_change_max = nb_change_max,distrib = distrib,mat_data = as.matrix(genomicmatrix),model = model)
    indices<-c(hicsegresults$t_hat)
    zerosremoved<-indices[!indices == 0]
    #truncating the block indicies for the shorter dimension
    #returning the row and column breakpoints (instead of just one set of them, as HiCseg typically does).
    return(zerosremoved)
  } else 
    {
    #extending the matrix for rows
    if(ncol(genomicmatrix) > nrow(genomicmatrix))
    {
      if(MI_strategy=="average") {    extended_matrix<-rbind(genomicmatrix,rep(rep(mean(as.numeric(unlist(genomicmatrix)))),ncol(genomicmatrix)),c(ncol(genomicmatrix)-nrow(genomicmatrix)))  }
      if(MI_strategy=="copy") {  extended_matrix<-rbind(genomicmatrix,genomicmatrix[1:(ncol(genomicmatrix)-nrow(genomicmatrix)),])}
      if(transpose){
        if(MI_strategy=="average") {    t_extended_matrix<-cbind(t(genomicmatrix),
                      matrix(rep(rep(mean(as.numeric(unlist(t(genomicmatrix)))),ncol(t(genomicmatrix))),c(nrow(t(genomicmatrix))-ncol(t(genomicmatrix)))),nrow = nrow(t(genomicmatrix)),ncol=c(nrow(t(genomicmatrix))-ncol(t(genomicmatrix)))
                      ))
        }
        if(MI_strategy=="copy") {  t_extended_matrix<-cbind(t(genomicmatrix),t(genomicmatrix)[,1:(nrow(t(genomicmatrix))-ncol(t(genomicmatrix)))]) 
        }
      }
      
    }
    if(ncol(genomicmatrix) < nrow(genomicmatrix))
    {
      if(MI_strategy=="average") {   extended_matrix<-cbind(genomicmatrix,rep(rep(mean(as.numeric(unlist(genomicmatrix)))),nrow(genomicmatrix)),c(nrow(genomicmatrix)-ncol(genomicmatrix)))  }
      if(MI_strategy=="copy") { extended_matrix<-cbind(genomicmatrix,genomicmatrix[,1:(nrow(genomicmatrix)-ncol(genomicmatrix))])  }
      if(transpose){
        if(MI_strategy=="average") {    t_extended_matrix<-rbind(t(genomicmatrix),matrix(rep(
          rep(mean(as.numeric(unlist(t(genomicmatrix)))) #a single mean
          ,nrow(t(genomicmatrix))) #a single row
          ,(ncol(t(genomicmatrix))-nrow(t(genomicmatrix)))),ncol=ncol(t(genomicmatrix)),nrow=(ncol(t(genomicmatrix))-nrow(t(genomicmatrix))))) }
          #the missing piece  
        if(MI_strategy=="copy") {  t_extended_matrix<-rbind(t(genomicmatrix),t(genomicmatrix)[1:(ncol(t(genomicmatrix))-nrow(t(genomicmatrix))),]) }
      }
    }
    if(transpose){
    t_hicsegresults<-HiCseg::HiCseg_linkC_R(size_mat=dim(t_extended_matrix)[1],nb_change_max = nb_change_max,distrib = distrib,mat_data = as.matrix(t_extended_matrix),model = model)
    t_indices<-c(t_hicsegresults$t_hat)
    t_zerosremoved<-t_indices[!t_indices == 0]
    t_breakpoints_col<-t_zerosremoved[t_zerosremoved < ncol(t(genomicmatrix))]
    t_breakpoints_row<-t_zerosremoved[t_zerosremoved < nrow(t(genomicmatrix))]
    }
      hicsegresults<-HiCseg::HiCseg_linkC_R(size_mat=dim(extended_matrix)[1],nb_change_max = nb_change_max,distrib = distrib,mat_data = as.matrix(extended_matrix),model = model)
  }
  indices<-c(hicsegresults$t_hat)
  zerosremoved<-indices[!indices == 0]
  #truncating the block indicies for the shorter dimension
  breakpoints_col<-zerosremoved[zerosremoved < ncol(genomicmatrix)]
  breakpoints_row<-zerosremoved[zerosremoved < nrow(genomicmatrix)]
  #returning the row and column breakpoints (instead of just one set of them, as HiCseg typically does).
  output_list<-list(breakpoints_col,breakpoints_row) #it's computationally more efficient to return these as a pair as they will almost always be needed in pairs.
  names(output_list)<-c("breakpoints_col","breakpoints_row")
  if(transpose) {output_list<-list(breakpoints_col,breakpoints_row,t_breakpoints_col,t_breakpoints_row)
  output_list<-list(breakpoints_col,breakpoints_row,t_breakpoints_col,t_breakpoints_row)
  names(output_list)<-c("breakpoints_col","breakpoints_row","t_breakpoints_col","t_breakpoints_row")
  }
  #genomicmatrix<-t(genomicmatrix)
  return(output_list)
}
