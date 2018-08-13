#'  Postprocess linear regression matrix.
#'
#' Takes a linear regression matrix and sets infinites to a finite value, and changes the sign to match the sign of the correlation for each value.
#' @keywords lm linear regression matrix correlation
#' @import matrixStats stats
#' @param input_matrix The input matrix, which consists of bins and samples (no LM or correlation has been done on the segmentation values)
#' @param LM_mat The linear regression matrix, with rows and columns consisting of bins and the values being the negative log p-value between them.
#' @param cor_type The correlation type ("pearson" (linear), "spearman" (rank), "kendall"(also rank-based)).
#'  Rank correlations capture nonlinear relationships as well as linear. Passed to stats::cor's method parameter.
#' @return The output matrix, or if using slurm, the slurm job object (which should be saved as an rds file and reloaded when creating the output matrix).
#' @examples
#' inputmat<-matrix(runif(15),nrow=3)
#' colnames(inputmat)<-c("chr2_1_1000","chr2_1001_2000","chr2_2001_3000","chr2_3001_4000","chr2_4001_5000")
#' rownames(inputmat)<-c("PAFPJK","PAKKAT","PUFFUM")
#' outputmat<-matrix(runif(15),nrow=3)
#' outputmat<-cor(inputmat)*matrix(runif(25,-30,500),nrow=5)
#' diag(outputmat)<-Inf
#' library(HiCNV)
#' postProcessLinRegMatrix(input_matrix=t(inputmat),LM_mat=outputmat,cor_type="pearson",inf_replacement_val=300)
#' @export
postProcessLinRegMatrix<-function(input_matrix,LM_mat,cor_type="pearson",inf_replacement_val=300)
{
  #removing empty columns:"
  #input_matrix_zeros_removed<-as.data.frame(t(input_matrix))[,colSums(as.data.frame(t(input_matrix)))>0]
  input_matrix_zeros_removed<-as.data.frame(t(input_matrix))[,colSds(as.matrix(t(input_matrix)))!=0] #this will take care of zero bins and invariant bins.
  #input_matrix_zeros<-as.data.frame(t(input_matrix))[,colSums(as.data.frame(t(input_matrix)))==0]
  #correcting infinites
  LM_mat[is.infinite(unlist(LM_mat))]<-inf_replacement_val
  TCGA_low_pass_matrix<-LM_mat
  #input_matrix<-input_matrix_zeros_removed
  #adding column names
  #this has to take account of the constant columns that are removed in the LM process 
  #These will not be represented as a linear regression p-value does not exist for two vectors where one is all the same value.
  rownames(TCGA_low_pass_matrix)<-colnames(input_matrix_zeros_removed)
  colnames(TCGA_low_pass_matrix)<-colnames(input_matrix_zeros_removed)
  #removes unnecessary substructure of the matrix/df.
  TCGA_low_pass_matrix2<-  matrix(as.numeric(unlist(TCGA_low_pass_matrix)),ncol=ncol(TCGA_low_pass_matrix))
  colnames(TCGA_low_pass_matrix2)<-colnames(TCGA_low_pass_matrix)
  rownames(TCGA_low_pass_matrix2)<-rownames(TCGA_low_pass_matrix)
  #fixing sign
  input_matrix_mat<-as.matrix(input_matrix_zeros_removed)
  input_matrix_cor<-cor(as.matrix(input_matrix_mat), use = "pairwise.complete.obs",method=cor_type)
  input_matrix_cor[input_matrix_cor>0]<-1
  input_matrix_cor[input_matrix_cor<0]<- -1
  TCGA_low_pass_matrix_sign_corrected<-(TCGA_low_pass_matrix2*input_matrix_cor)
  return(TCGA_low_pass_matrix_sign_corrected)
}