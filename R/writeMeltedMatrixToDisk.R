#' Write a matrix, with genes, of a submatrix of a whole genome interaction matrix to disk.
#'
#' Writes an RData file with a ggplot2 object within.
#' @keywords ggplot2 heatmaply plotly ggiraph genomic matrix
#' @import IRanges GenomicRanges magrittr OpenImageR
#' @param whole_genome_matrix A matrix to have edges averaged with genomic coordinates in the form chr1_50_100 set as the column and row names.
#' @param chrom1 first chromosome of the two which will subset the matrix. (this is done in row-column fasion).
#' @param chrom2 second chromosome of the two which will subset the matrix. (this is done in row-column fasion).
#' @param extra_data_matrix A matrix with additional variables about each point, one position per row with as many variables as remaining columns.
#' @return ggplotmatrix a matrix with values sufficient to create a ggplot2 heatmap with geom_tile() or with ggiraph's geom_tile_interactive()
#' @export
writeMeltedChromosomalMatrixToDisk<-function(whole_genome_matrix,chrom1,chrom2,filename,extra_data_matrix=NULL,transpose=F,sequential=T,debug=T,multipass=T,desired_range_start=50,desired_range_end=300)
{
  if(!is.null(extra_data_matrix))
  {  
    extra_data_matrix_fn<-extra_data_matrix
    extra_data_df<-data.table::fread(extra_data_matrix,data.table=F)
  }
  chromosomes<-paste0("chr",c(seq(1:22),"X"),"_")
  submatrix<-whole_genome_matrix[grep(chromosomes[chrom1],rownames(whole_genome_matrix)),grep(chromosomes[chrom2],colnames(whole_genome_matrix))]
  #insert intra code here for compatibility, remembering the bit at the end of the while loops.
  downsample_factor<-NULL
  desired_range<-IRanges(desired_range_start,desired_range_end)
  downsample_factor_row<-NULL
  downsample_outcomes_row<-as.data.frame(cbind(numbers::divisors(nrow(submatrix)),nrow(submatrix)/numbers::divisors(nrow(submatrix))))
  colnames(downsample_outcomes_row)<-c("factor","downsampled_size")
  downsample_factor_col<-NULL
  downsample_outcomes_col<-as.data.frame(cbind(numbers::divisors(ncol(submatrix)),ncol(submatrix)/numbers::divisors(ncol(submatrix))))
  colnames(downsample_outcomes_col)<-c("factor","downsampled_size")
  while(length(intersect(downsample_factor_col,downsample_factor_row))==0 & (nrow(submatrix)>desired_range_start & ncol(submatrix)>desired_range_start))
  {
    downsample_factor<-NULL
    downsample_factor_row<-NULL
    downsample_factor_col<-NULL
    downsample_outcomes<-NULL
    while(length(downsample_factor_col)==0)
    {
      downsample_outcomes_row<-as.data.frame(cbind(numbers::divisors(nrow(submatrix)),nrow(submatrix)/numbers::divisors(nrow(submatrix))))
      colnames(downsample_outcomes_row)<-c("factor","downsampled_size")
      downsample_factor_row<-downsample_outcomes_row[downsample_outcomes_row$downsampled_size>=desired_range@start & downsample_outcomes_row$downsampled_size<=(desired_range@start+desired_range@width),"factor"]
      downsample_outcomes_col<-as.data.frame(cbind(numbers::divisors(ncol(submatrix)),ncol(submatrix)/numbers::divisors(ncol(submatrix))))
      colnames(downsample_outcomes_col)<-c("factor","downsampled_size")
      downsample_factor_col<-downsample_outcomes_col[downsample_outcomes_col$downsampled_size>=desired_range@start & downsample_outcomes_col$downsampled_size<=(desired_range@start+desired_range@width),"factor"]
      submatrix_temp<-submatrix
      if(debug){    print(paste0("col factors:",downsample_factor_col))
        print(paste0("row factors in desired range:",downsample_factor_row))
        print(paste0("current submatrix dimensions:",paste0(dim(submatrix))))
        print(paste0("row factors (including those outside desired range):",paste(downsample_outcomes_row$factor)))
        print(paste0("col factors (including those outside desired range):",paste(downsample_outcomes_col$factor)))}
      #browser()
      #temp_ds_factor<-as.integer(grep(1,intersect(downsample_outcomes_row$factor,downsample_outcomes_col$factor),value=T,invert=T))[1]
      if(length(intersect(downsample_outcomes_row$factor,downsample_outcomes_col$factor))>1)  {submatrix<-downsample_genomic_matrix(whole_matrix=submatrix,downsamplefactor=as.integer(grep(1,intersect(downsample_outcomes_row$factor,downsample_outcomes_col$factor),value=T,invert=T))[1])} else {
        
      # if(length(intersect(downsample_outcomes_row$factor,downsample_outcomes_col$factor))>1)  {
      #   temp_ds_factor<-as.integer(grep(1,intersect(downsample_outcomes_row$factor,downsample_outcomes_col$factor),value=T,invert=T))[1]
      #   if(ncol(submatrix)/temp_ds_factor>=desired_range_start & nrow(submatrix)/temp_ds_factor>=desired_range_start){
      #   submatrix<-downsample_genomic_matrix(whole_matrix=submatrix,downsamplefactor=temp_ds_factor)
      #   }
      #   rm(temp_ds_factor)
      #   } else {
      if(length(downsample_factor_col)==0 & ncol(submatrix)>(desired_range@start+1)){submatrix<-averageMatrixEdges(submatrix,dimension = "column")}
      } #only average edges if there isn't another way to downsample first. #I should change this such that it doens't downsample too much.  This will likely be another edge case.
    }
    
    #downsample_outcomes_col<-downsample_outcomes
    #downsample_factor_col<-downsample_factor
    downsample_factor<-NULL
if(debug){    print(paste0("col factors:",downsample_factor_col))
    print(paste0("row factors:",downsample_factor_row))
    print(paste0("current submatrix dimensions:",paste0(dim(submatrix)))) }
    while(length(downsample_factor_row)==0)
    {
      downsample_outcomes_col<-as.data.frame(cbind(numbers::divisors(ncol(submatrix)),ncol(submatrix)/numbers::divisors(ncol(submatrix))))
      colnames(downsample_outcomes_col)<-c("factor","downsampled_size")
      downsample_factor_col<-downsample_outcomes_col[downsample_outcomes_col$downsampled_size>=desired_range@start & downsample_outcomes_col$downsampled_size<=(desired_range@start+desired_range@width),"factor"]
      downsample_outcomes_row<-as.data.frame(cbind(numbers::divisors(nrow(submatrix)),nrow(submatrix)/numbers::divisors(nrow(submatrix))))
      colnames(downsample_outcomes_row)<-c("factor","downsampled_size")
      downsample_factor_row<-downsample_outcomes_row[downsample_outcomes_row$downsampled_size>=desired_range@start & downsample_outcomes_row$downsampled_size<=(desired_range@start+desired_range@width),"factor"]
      submatrix_temp<-submatrix
      if(debug){    print(paste0("row factors:",downsample_factor_row))
        print(paste0("col factors:",downsample_factor_col))
        print(paste0("row factors in desired range:",downsample_factor_row))
        print(paste0("current submatrix dimensions:",paste0(dim(submatrix))))
        print(paste0("row factors (including those outside desired range):",paste(downsample_outcomes_row$factor)))
        print(paste0("col factors (including those outside desired range):",paste(downsample_outcomes_col$factor)))}
      #browser()
           if(length(intersect(downsample_outcomes_row$factor,downsample_outcomes_col$factor))>1)  {submatrix<-downsample_genomic_matrix(whole_matrix=submatrix,downsamplefactor=as.integer(grep(1,intersect(downsample_outcomes_row$factor,downsample_outcomes_col$factor),value=T,invert=T))[1])} else {

      # if(length(intersect(downsample_outcomes_row$factor,downsample_outcomes_col$factor))>1)  {
      #   temp_ds_factor<-as.integer(grep(1,intersect(downsample_outcomes_row$factor,downsample_outcomes_col$factor),value=T,invert=T))[1]
      #   if(ncol(submatrix)/temp_ds_factor>=desired_range_start & nrow(submatrix)/temp_ds_factor>=desired_range_start){
      #     submatrix<-downsample_genomic_matrix(whole_matrix=submatrix,downsamplefactor=temp_ds_factor)
      #   }
      #   rm(temp_ds_factor)
      # } else {        
      if(length(downsample_factor_row)==0 & nrow(submatrix)>(desired_range@start+1)){submatrix<-averageMatrixEdges(submatrix,dimension = "row")}
      }
    }
    #downsample_outcomes_row<-downsample_outcomes
    #downsample_factor_row<-downsample_factor
    #browser()
    downsample_outcomes_col<-as.data.frame(cbind(numbers::divisors(ncol(submatrix)),ncol(submatrix)/numbers::divisors(ncol(submatrix))))
    colnames(downsample_outcomes_col)<-c("factor","downsampled_size")
    downsample_factor_col<-downsample_outcomes_col[downsample_outcomes_col$downsampled_size>=desired_range@start & downsample_outcomes_col$downsampled_size<=(desired_range@start+desired_range@width),"factor"]
    downsample_outcomes_row<-as.data.frame(cbind(numbers::divisors(nrow(submatrix)),nrow(submatrix)/numbers::divisors(nrow(submatrix))))
    colnames(downsample_outcomes_row)<-c("factor","downsampled_size")
    downsample_factor_row<-downsample_outcomes_row[downsample_outcomes_row$downsampled_size>=desired_range@start & downsample_outcomes_row$downsampled_size<=(desired_range@start+desired_range@width),"factor"]
    submatrix_temp<-submatrix
    #browser()
    # if(length(intersect(downsample_outcomes_row$factor,downsample_outcomes_col$factor))>1)  {
    #   temp_ds_factor<-as.integer(grep(1,intersect(downsample_outcomes_row$factor,downsample_outcomes_col$factor),value=T,invert=T))[1]
    #   if(ncol(submatrix)/temp_ds_factor>=desired_range_start & nrow(submatrix)/temp_ds_factor>=desired_range_start){
    #     submatrix<-downsample_genomic_matrix(whole_matrix=submatrix,downsamplefactor=temp_ds_factor)
    #   }
    #   rm(temp_ds_factor)
    # } else {
    if(length(intersect(downsample_outcomes_row$factor,downsample_outcomes_col$factor))>1)  {submatrix<-downsample_genomic_matrix(whole_matrix=submatrix,downsamplefactor=as.integer(grep(1,intersect(downsample_outcomes_row$factor,downsample_outcomes_col$factor),value=T,invert=T))[1])} else {
      
#     if(length(intersect(downsample_outcomes_row$factor,downsample_outcomes_col$factor))>1)  {submatrix<-downsample_genomic_matrix(whole_matrix=submatrix,downsamplefactor=as.integer(grep(1,intersect(downsample_outcomes_row$factor,downsample_outcomes_col$factor),value=T,invert=T))[1])} else {
    if(length(intersect(downsample_factor_col,downsample_factor_row))==0) {
      if(length(downsample_factor_row)>length(downsample_factor_col)){ submatrix<-averageMatrixEdges(submatrix,dimension="row")} else {
        if(length(downsample_factor_row)<length(downsample_factor_col)){submatrix<-averageMatrixEdges(submatrix,dimension="column")}
        if(length(downsample_factor_row)==length(downsample_factor_col) & nrow(downsample_outcomes_row)>=nrow(downsample_outcomes_col)){submatrix<-averageMatrixEdges(submatrix,dimension="column")}
        if(length(downsample_factor_row)==length(downsample_factor_col) & nrow(downsample_outcomes_row)<nrow(downsample_outcomes_col)){submatrix<-averageMatrixEdges(submatrix,dimension="column")}
        }#end else
    }
    } #end outer else
    if(debug){print(paste0("col factors:",downsample_factor_col))
    print(paste0("row factors:",downsample_factor_row))
    print(paste0("current submatrix dimensions:",paste0(dim(submatrix))))}
  } #end while
  if(nrow(submatrix)==ncol(submatrix)){if(length(downsample_factor_col)>0){downsample_factor<-min(downsample_factor_col)}}
  if(nrow(submatrix)!=ncol(submatrix)){if(length(intersect(downsample_factor_row,downsample_factor_col))>0){downsample_factor<-min(intersect(downsample_factor_row,downsample_factor_col))}}
 if(debug){ print(paste0("final col:",downsample_factor_col)) 
   print(paste0("final row:",downsample_factor_row))}
  
  submatrix_downsample<-downsample_genomic_matrix(submatrix,downsample_factor,singlechromosome = T)
  if(transpose){concatenated_gene_matrix<-getAnnotationMatrix(t(submatrix_downsample),prot_only = T,flip_row_col=T,sequential=T)} else{
    concatenated_gene_matrix<-getAnnotationMatrix(submatrix_downsample,prot_only = T,flip_row_col=T)}
  concatenated_gene_matrix.m<-melt(concatenated_gene_matrix)
  concatenated_gene_matrix.m$Var1<-rownames(submatrix_downsample)[concatenated_gene_matrix.m$Var1]
  concatenated_gene_matrix.m$Var2<-colnames(submatrix_downsample)[concatenated_gene_matrix.m$Var2]
  if(transpose){
    ggplotmatrix<-t(submatrix_downsample) %>% as.matrix() %>%  signedRescale() %>%
      melt() %>% dplyr::bind_cols(concatenated_gene_matrix.m)
    
  } else {
    ggplotmatrix<-submatrix_downsample %>% as.matrix() %>%  signedRescale() %>%
      melt() %>% dplyr::bind_cols(concatenated_gene_matrix.m)
  }
  if(!is.null(extra_data_matrix))
  {
    #need to also downsample the extra_data_df to match the positions of the ggplotmatrix$Var1 and Var2 (which should be equal).
    ggplotmatrix$Var1<-as.character(ggplotmatrix$Var1)
    ggplotmatrix$Var2<-as.character(ggplotmatrix$Var2)
    row_merged_ggplotmatrix<-merge(ggplotmatrix,extra_data_df,by.x="Var1",by.y="pos",suffixes=c("","row")) #var1 row Var2 Column, y x.
    row_col_merged_ggplotmatrix<-merge(row_merged_ggplotmatrix,extra_data_df,by.x="Var2",by.y="pos",suffixes=c("","col"))
    print(head(row_col_merged_ggplotmatrix,n=1))
    if(transpose){save("row_col_merged_ggplotmatrix",file=paste0(chromosomes[chrom1],chromosomes[chrom2],"melted_with_external_data_transposed",".RData"))} else {save("row_col_merged_ggplotmatrix",file=paste0(chromosomes[chrom1],chromosomes[chrom2],"melted_with_external_data",".RData"))}
    
    
  } else {
    #write.csv(ggplotmatrix,paste0(chromosomes[chrom1],chromosomes[chrom2],"melted",".csv"),row.names = F)
    save("ggplotmatrix",file=paste0(chromosomes[chrom1],chromosomes[chrom2],"melted",".RData"))
    #rm("ggplotmatrix")
  }
  return(ggplotmatrix)
}
