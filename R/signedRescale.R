#' Rescale positive and negative data, preserving sign information.
#'
#' Performs a signed rescale on the data, shrinking the negative and positive ranges into the [0,1] space, such that negative is always less than 0.5 and positive is always greater.
#' @keywords signed rescale positive negative matrix
#' @param matrix A matrix to be transformed
#' @param global_max the global maximum (used if scaling using statistics from a large matrix upon a submatrix).
#' @param global_min the global minimum
#' @param global_sigma the global signma
#' @param global_mu the global mu
#' @param max_cap the maximum saturation-- decreases the ceiling considered for the scaling function.
#' Useful to see greater differences if an image is too white, increase it if there is too much color to tell apart domains.
#' @param method method to perform the rescaling.
#' Options are "minmax" (default), "tan" for tangent, and "sd" for standard devation
#' @param tan_transform apply a tangent transformation?
#' @param global_sigma_pos The positive global sigma. See getGlobalRescalingStats. 
#' @param global_sigma_neg The negative global sigma. See getGlobalRescalingStats.
#' @param asymptotic_max make the maximum value in the matrix not 1, but rather something slightly below.
#' @return transformedmatrix A transformed matrix.
#' @examples 
#' \dontrun{
#' mat<-matrix(c(5,10,15,20,0,40,-45,300,-50),byrow=TRUE,nrow=3)
#' rescaled_mat<-signedRescale(mat)
#' mat
#' rescaled_mat
#' }
#' @export
signedRescale<-function(matrix,global_max=NULL,global_min=NULL,global_sigma=NULL,global_mu=NULL,max_cap=NULL,method="minmax",tan_transform=F,global_sigma_pos=NULL,global_sigma_neg=NULL,asymptotic_max=T)
{
  #matrix<-as.matrix(matrix)
  transformedmatrix<-as.matrix(matrix)
  transformedmatrix[transformedmatrix==0]<-(0.5+1e9*.Machine$double.eps)
  if(!is.null(max_cap)){transformedmatrix[transformedmatrix>max_cap]<-max_cap}
  if(is.null(global_max)){global_max<-max(transformedmatrix[transformedmatrix>0])}
  if(is.null(global_min)){global_min<-min(transformedmatrix[transformedmatrix<0])}
  if(is.null(global_sigma_pos)){global_sigma_pos<-sd(transformedmatrix[transformedmatrix>0])}
  if(is.null(global_sigma_neg)){global_sigma_neg<-sd(transformedmatrix[transformedmatrix<0])}
  if(is.null(global_sigma)){global_sigma<-sd(transformedmatrix)}
  if(is.null(global_mu)){global_mu<-mean(transformedmatrix)}
  
  
  #if(tan_transform){transformedmatrix<-atan(transformedmatrix)/pi*2}
  if(method=="minmax"){
    #browser()
    transformedmatrix[transformedmatrix>0 & transformedmatrix!=(0.5+1e9*.Machine$double.eps)]<-((transformedmatrix[transformedmatrix>0 & transformedmatrix!=(0.5+1e9*.Machine$double.eps)]/(global_max*2))+0.5) #divide by global max * 2, store into transformed matrix.
    transformedmatrix[matrix<0]<-((transformedmatrix[transformedmatrix<0]/(global_min*2))) #divide by global minimum * 2, store into transformed matrix.
    transformedmatrix[transformedmatrix<=0.5 & transformedmatrix>0]<-abs(0.5-transformedmatrix[transformedmatrix<=0.5 & transformedmatrix>0]) #(0,0.5), negative numbers in original matrix.
    #transformedmatrix[transformedmatrix==0]<-(0.5-1e9*.Machine$double.eps)
  }
  if(method=="tan"){transformedmatrix<-atan(transformedmatrix)/pi*2}
  if(method=="sd")
  {
    transformedmatrix[transformedmatrix>0]<-((transformedmatrix[transformedmatrix>0]/(global_max*2))+0.5)
    transformedmatrix[matrix<0]<-((transformedmatrix[transformedmatrix<0]/(global_min*2)))
    transformedmatrix[transformedmatrix<0.5 & transformedmatrix>0]<-abs(0.5-transformedmatrix[transformedmatrix<0.5 & transformedmatrix>0])
    transformedmatrix[transformedmatrix==0]<-(0.5-1e9*.Machine$double.eps)
  }
  if(asymptotic_max)
  {
    #browser()
    transformedmatrix[transformedmatrix==1]<-(1-1e9*.Machine$double.eps)
  }
  return(transformedmatrix)
}
#all_conc_cleaned_common_coords_linreg_tiny<-all_conc_cleaned_common_coords_linreg[1:25,1:25]
# all_conc_cleaned_common_coords_linreg_tiny.m<-melt(as.matrix(all_conc_cleaned_common_coords_linreg[1:5,1:5]))
# signedRescale(all_conc_cleaned_common_coords_linreg_tiny)==HiCNV::signedRescale(all_conc_cleaned_common_coords_linreg_tiny)
# bins<-fread(paste0(groupdir,"dalgleishjl/hicnv/binfile.txt"))
# bins_t<-t(bins[,2:ncol(bins)])
# colnames(bins_t)<-bins$probe
# #head(bins_t)
# Heatmap(signedRescale(all_conc_cleaned_common_coords_linreg_tiny,tan_transform = F),cluster_rows = F,cluster_columns = F,colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")))
# Heatmap(HiCNV::signedRescale(as.matrix(all_conc_cleaned_common_coords_linreg[1:1107,1:1107]),tan_transform = F),cluster_rows = F,cluster_columns = F,colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")))
# Heatmap(signedRescale(as.matrix(all_conc_cleaned_common_coords_linreg[1:1107,1:1107]),tan_transform = F),cluster_rows = F,show_row_names=F,show_column_names=F,cluster_columns = F,colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")))
# HiCNV::signedRescale(as.matrix(all_conc_cleaned_common_coords_linreg[1:85,1:85]))
# matrix=all_conc_cleaned_common_coords_linreg[1:85,1:85]
# setwd(paste0(groupdir,"dalgleishjl/hicnv/color_scale_test_plots/"))
# 
# Heatmap(signedRescale(as.matrix(all_conc_cleaned_common_coords_linreg[1:1107,1:1107]),method="tan",tan_transform = F),cluster_rows = F,show_row_names=F,show_column_names=F,cluster_columns = F,colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")))
# 
# Heatmap(signedRescale(as.matrix(all_conc_cleaned_common_coords_linreg[1:1107,1:1107]),method="minmax",tan_transform = F,max_cap = 200),cluster_rows = F,show_row_names=F,show_column_names=F,cluster_columns = F,colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")))
# Heatmap(signedRescale(as.matrix(all_conc_cleaned_common_coords_linreg[1:1107,1:1107]),method="minmax",tan_transform = F,max_cap = 100),cluster_rows = F,show_row_names=F,show_column_names=F,cluster_columns = F,colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")))
# Heatmap(signedRescale(as.matrix(all_conc_cleaned_common_coords_linreg[1:1107,1:1107]),method="minmax",tan_transform = F,max_cap = 50),cluster_rows = F,show_row_names=F,show_column_names=F,cluster_columns = F,colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")))
# Heatmap(signedRescale(as.matrix(all_conc_cleaned_common_coords_linreg[1:1107,1:1107]),method="minmax",tan_transform = F,max_cap = 25),cluster_rows = F,show_row_names=F,show_column_names=F,cluster_columns = F,colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")))
# foreach(i=c(10,25,50,75,100,200)) %dopar%
# {
#   png(paste0("chr1_max_cap",i,".png"))
#   print(Heatmap(signedRescale(as.matrix(all_conc_cleaned_common_coords_linreg[1:1107,1:1107]),method="minmax",tan_transform = F,max_cap = i),cluster_rows = F,show_row_names=F,show_column_names=F,cluster_columns = F,colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))))
#   dev.off()
# }
# png(paste0("chr1_tan_transform",".png"))
# print(Heatmap(signedRescale(as.matrix(all_conc_cleaned_common_coords_linreg[1:1107,1:1107]),method="tan",tan_transform = F),cluster_rows = F,show_row_names=F,show_column_names=F,cluster_columns = F,colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))))
# dev.off()