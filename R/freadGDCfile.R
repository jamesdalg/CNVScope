freadGDCfile<-function(file,fread_skip=NULL) {input_csv<-fread(file,skip=fread_skip)
sample_info_colsplit<-reshape2::colsplit(basename(file),"_|-|\\.",c("pre","project","num","sample","comparison","fn_ext"))
input_csv_with_sample_info<-dplyr::bind_cols(input_csv,sample_info_colsplit[rep(1,nrow(input_csv)),])
if(length(na.omit(unlist(sample_info_colsplit)))!=6){return(NULL)}
return(input_csv_with_sample_info)
}
