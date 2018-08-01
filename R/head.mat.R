#' Gets a small piece of a matrix (top left corner) for viewing, rather than pulling the first n rows.
#'
#' Gives a small square of a matrix to get an idea of content rather than grabbing the entire row.
#' When this row is thousands of numbers long, this can be a problem.
#' @keywords rescale downsample average edges matrix
#' @param mat A matrix.
#' @param n The length and width of the piece to view.
#' @return averaged_matrix a small matrix of size n.
#' @export

head.mat<-function(mat,n=6L)
{
  return(mat[1:n,1:n])
}
