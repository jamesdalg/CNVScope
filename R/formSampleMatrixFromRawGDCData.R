#' Form sample matrix from GDC low-pass segmentation datafiles.
#'
#' Reads a GDC segmetnation files, adds sample information, and forms a data matrix of samples and bins of a specified size.
#' @keywords segmentation GDC 
#' @import reshape2 dplyr data.table BSgenome.Hsapiens.UCSC.hg19 doMC
#' @param file GDC file to be read
#' @param format file format, TCGA or TARGET.
#' @param binsize the binsize, in base pairs (default 1Mb or 1e6).  This value provides a good balance of resolution and speed with memory sensitive applications.
#' @param freadksip the number of lines to skip in the GDC files, typically 14 (the first 13 lines are metadata and the first is a blank line in NBL data). Adjust as needed.
#' @return sample_aggregated_segvals A dataframe containing the aggregated segmentation values, based on the parameters provided.
#' @export


formSampleMatrixFromRawGDCData<-function(tcga_files=NULL,format="TARGET",binsize=1e6,freadskip=14)
{
  chromosomes<-paste0("chr",c(seq(1:22),"X"),"_")
  # TCGA_CNV_data_with_sample_info<-ldply(tcga_files,
  #                                       function(x) {input_csv<-fread(x,skip=freadskip)
  #                                       sample_info_colsplit<-reshape2::colsplit(basename(x),"_|-|\\.",c("pre","project","num","sample","comparison","fn_ext"))
  #                                       input_csv_with_sample_info<-dplyr::bind_cols(input_csv,sample_info_colsplit[rep(1,nrow(input_csv)),])
  #                                       return(input_csv_with_sample_info)
  #                                       }
  # )
  TCGA_CNV_data_with_sample_info<-ldply(tcga_files,function(x) freadGDCfile(x,fread_skip=freadskip))  
  #TCGA_CNV_data_with_sample_info_small<-ldply(tcga_files[1:100],freadGDCfile)  
  #end testing
  TCGA_CNV_data<-TCGA_CNV_data_with_sample_info
  if(format=="TCGA"){
    TCGA_CNV_data$Chromosome<-paste0("chr",TCGA_CNV_data$Chromosome)
    colnames(TCGA_CNV_data)<-gsub("Chromosome",">chr",gsub("End","end",gsub("Start","begin",colnames(TCGA_CNV_data))))
  }
  TCGA_CNV_data_range_filtered<-TCGA_CNV_data %>% tidyr::drop_na(begin,end)
  TCGA_CNV_data_dt<-as.data.table(TCGA_CNV_data_range_filtered)
  TCGA_CNV_data_gr<-GRanges(seqnames = TCGA_CNV_data_range_filtered$`>chr`,ranges = IRanges(start = TCGA_CNV_data_range_filtered$begin,end = TCGA_CNV_data_range_filtered$end),... = TCGA_CNV_data_range_filtered[,4:ncol(TCGA_CNV_data_range_filtered)])
  bins<-tileGenome(seqinfo(Hsapiens),tilewidth=binsize,cut.last.tile.in.chrom = T)
  bins<-bins[bins@seqnames %in% gsub("_","",chromosomes)]
  rownames_gr = bins
  colnames_gr = bins
  samples<-unique(mcols(TCGA_CNV_data_gr)$....sample)
  options(scipen=999)
  bins_underscored<-GRanges_to_underscored_pos(bins)
  registerDoMC()
  TCGA_CNV_data_gr_all_comparisons<-TCGA_CNV_data_gr
  TCGA_CNV_data_gr_single_comparison<-TCGA_CNV_data_gr[mcols(TCGA_CNV_data_gr)$....comparison=="NormalVsPrimary"]
  TCGA_CNV_data_gr<-TCGA_CNV_data_gr_single_comparison
  sample_aggregated_segvals<-foreach(s=1:length(samples),.combine="cbind",.errorhandling = "stop") %dopar% {
    current_gr<-TCGA_CNV_data_gr[mcols(TCGA_CNV_data_gr)$....sample %in% samples[s]]
    current_merged_df<-as.data.frame(mergeByOverlaps(bins,current_gr))
    current_merged_df$pos<-unlist(tidyr::unite(current_merged_df[,c("bins.seqnames","bins.start","bins.end")]))
    #sort(table(current_merged_df$pos),decreasing=T)
    current_merged_df_bins_vals<-current_merged_df[,c("pos","....relativeCvg","....sample")] #,"....comparison"
    current_merged_df_bins_vals$....relativeCvg<-as.numeric(as.character(current_merged_df_bins_vals$....relativeCvg))
    current_merged_df_bins_vals<-na.omit(current_merged_df_bins_vals)
    #current_merged_df_bins_aggregated_test<-ddply(na.omit(current_merged_df_bins_vals[1,]),.(pos),summarise,meanrelcvg=mean(current_merged_df_bins_vals$....relativeCvg))#,samples=list(unique(current_merged_df_bins_vals$....sample))
    
    current_merged_df_bins_aggregated<-ddply(na.omit(current_merged_df_bins_vals),.(pos),summarise,meanrelcvg=mean(....relativeCvg),samples=paste0(unique(....sample),collapse=","))#
    #insert bins that are not represented.
    unused_bins<-bins_underscored[!(bins_underscored %in% current_merged_df_bins_aggregated$pos)]
    unused_bins_rows<-as.data.frame(cbind(unused_bins,rep(0,length(unused_bins)),rep(samples[s],length(unused_bins))))
    colnames(unused_bins_rows)<-c("pos","meanrelcvg","samples")
    unused_bins_rows$meanrelcvg<-as.numeric(unused_bins_rows$meanrelcvg)
    current_merged_df_bins_aggregated_with_unused<-rbind(current_merged_df_bins_aggregated[,c("pos","meanrelcvg","samples")],unused_bins_rows[,c("pos","meanrelcvg","samples")])
    current_merged_df_bins_aggregated_with_unused<-current_merged_df_bins_aggregated_with_unused[order(underscored_pos_to_GRanges(current_merged_df_bins_aggregated_with_unused$pos)),]
    #end testing
    relcvg_df<-as.data.frame(current_merged_df_bins_aggregated_with_unused$meanrelcvg)
    rownames(relcvg_df)<-current_merged_df_bins_aggregated_with_unused$pos
    colnames(relcvg_df)<-samples[s]
    
    print(paste0(samples[s]," complete"))
    return(relcvg_df)
  }
  return(sample_aggregated_segvals)
}
