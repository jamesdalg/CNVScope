#' Form sample matrix from GDC low-pass segmentation datafiles.
#'
#' Reads a GDC segmetnation files, adds sample information, and forms a data matrix of samples and bins of a specified size.
#' @name formSampleMatrixFromRawGDCData
#' @keywords segmentation GDC 
#' @import doParallel
#' @importFrom data.table fread
#' @importFrom reshape2 colsplit
#' @importFrom tidyr drop_na unite
#' @importFrom GenomicRanges tileGenome mcols
#' @importFrom IRanges mergeByOverlaps IRanges
#' @importFrom plyr ddply
#' @param tcga_files GDC files to be read
#' @param format file format, TCGA or TARGET.
#' @param binsize the binsize, in base pairs (default 1Mb or 1e6).  This value provides a good balance of resolution and speed with memory sensitive applications.
#' @param freadskip the number of lines to skip in the GDC files, typically 14 (the first 13 lines are metadata and the first is a blank line in NBL data). Adjust as needed.
#' @return sample_aggregated_segvals A dataframe containing the aggregated segmentation values, based on the parameters provided.
#' @examples 
#' #Pipeline examples would be too large to include in package checks.
#' #please see browseVignettes("CNVScope") for a demonstration.
#' 
#' @export
globalVariables(c('begin','s',".","pos",'....relativeCvg','....sample'))

formSampleMatrixFromRawGDCData<-function(tcga_files=NULL,format="TARGET",binsize=1e6,freadskip=NULL, parallel = F)
{
  chromosomes<-paste0("chr",c(seq(1:22),"X"),"_")
  # TCGA_CNV_data_with_sample_info<-ldply(tcga_files,
  #                                       function(x) {input_csv<-fread(x,skip=freadskip)
  #                                       sample_info_colsplit<-reshape2::colsplit(basename(x),"_|-|\\.",c("pre","project","num","sample","comparison","fn_ext"))
  #                                       input_csv_with_sample_info<-dplyr::bind_cols(input_csv,sample_info_colsplit[rep(1,nrow(input_csv)),])
  #                                       return(input_csv_with_sample_info)
  #                                       }
  # )
  if(format=="TARGET")
  {
    if(is.null(freadskip)) {freadskip=14}
    TCGA_CNV_data_with_sample_info<-plyr::ldply(tcga_files,function(x) freadGDCfile(x,fread_skip=freadskip)) 
    #TCGA_CNV_data_with_sample_info <- tcga_files %>% purrr::map_dfr(freadGDCfile, format = "TCGA",fread_skip=freadskip)
  }
  
  if(format=="TCGA"){
    TCGA_CNV_data_with_sample_info<-plyr::ldply(tcga_files,freadGDCfile,format="TCGA") #TCGA format files are really simple to read, compared to TARGET files.
    #TCGA_CNV_data_with_sample_info$....sample<-tcga_files
    #TCGA_CNV_data_with_sample_info <- tcga_files %>% purrr::map_dfr(freadGDCfile, format = "TCGA",fread_skip=freadskip)
    #return(TCGA_CNV_data_with_sample_info)
  }
  
  #TCGA_CNV_data_with_sample_info_small<-ldply(tcga_files[1:100],freadGDCfile)  
  #end testing
  TCGA_CNV_data<-TCGA_CNV_data_with_sample_info
  if(format=="TCGA"){
    if(!(all(stringr::str_detect(TCGA_CNV_data$Chromosome,"chr")))){
      TCGA_CNV_data$Chromosome<-paste0("chr",TCGA_CNV_data$Chromosome)
    }
    colnames(TCGA_CNV_data)<-gsub("Chromosome",">chr",gsub("End","end",gsub("Start","begin",colnames(TCGA_CNV_data))))
  }
  
  TCGA_CNV_data_range_filtered<-TCGA_CNV_data %>% tidyr::drop_na(begin,end)
  TCGA_CNV_data_dt<-data.table::as.data.table(TCGA_CNV_data_range_filtered)
  TCGA_CNV_data_gr<-GenomicRanges::GRanges(seqnames = TCGA_CNV_data_range_filtered$`>chr`,ranges = IRanges::IRanges(start = TCGA_CNV_data_range_filtered$begin,end = TCGA_CNV_data_range_filtered$end),... = TCGA_CNV_data_range_filtered[,4:ncol(TCGA_CNV_data_range_filtered)])
  bins<-GenomicRanges::tileGenome(seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens),tilewidth=binsize,cut.last.tile.in.chrom = T)
  bins<-bins[bins@seqnames %in% gsub("_","",chromosomes)]
  rownames_gr = bins
  colnames_gr = bins
  if(format=="TARGET")
  {samples<-unique(GenomicRanges::mcols(TCGA_CNV_data_gr)$....sample)}
  if(format=="TCGA")
  {
    #samples <- dirname(tcga_files) %>% tibble::as.tibble() %>% tidyr::separate(value, sep="/",into=c("dir","uuid")) %>% dplyr::pull(uuid)
    #mcols(TCGA_CNV_data_gr)$sample<-samples #best not done here.
    samples<-unique(GenomicRanges::mcols(TCGA_CNV_data_gr)$....uuid)
    GenomicRanges::mcols(TCGA_CNV_data_gr)$....sample<-GenomicRanges::mcols(TCGA_CNV_data_gr)$....uuid
    #add in sample column, add in sample mcol
    }
  options(scipen=999)
  bins_underscored<-GRanges_to_underscored_pos(bins)
  if(parallel){doParallel::registerDoParallel()}
  if(!parallel){foreach::registerDoSEQ()}
  if(format=="TARGET"){
  TCGA_CNV_data_gr_all_comparisons<-TCGA_CNV_data_gr
  TCGA_CNV_data_gr_single_comparison<-TCGA_CNV_data_gr[mcols(TCGA_CNV_data_gr)$....comparison=="NormalVsPrimary"]
  
  TCGA_CNV_data_gr<-TCGA_CNV_data_gr_single_comparison
  }
  sample_aggregated_segvals<-foreach(s=1:length(samples),.combine="cbind",.errorhandling = "stop",.export=ls(),.packages=c("magrittr","GenomicRanges","plyr","CNVScope")) %dopar% {
    #browser()
    if(format=="TCGA"){
      current_gr<-TCGA_CNV_data_gr[mcols(TCGA_CNV_data_gr)$....uuid %in% samples[s]]
    } 
    if(format == "TARGET")
    {
      current_gr<-TCGA_CNV_data_gr[mcols(TCGA_CNV_data_gr)$....sample %in% samples[s]]
    }
    current_merged_df<-as.data.frame(IRanges::mergeByOverlaps(bins,current_gr))
    current_merged_df$pos<-unlist(
      tidyr::unite(
      current_merged_df[,c("bins.seqnames","bins.start","bins.end")],pos
      )
      )
    #sort(table(current_merged_df$pos),decreasing=T)
    if(format=="TARGET")
    {
    current_merged_df_bins_vals<-current_merged_df[,c("pos","....relativeCvg","....sample")] #,"....comparison"
    current_merged_df_bins_vals$....relativeCvg<-as.numeric(as.character(current_merged_df_bins_vals$....relativeCvg))
    current_merged_df_bins_vals<-na.omit(current_merged_df_bins_vals)
    current_merged_df_bins_aggregated<-plyr::ddply(na.omit(current_merged_df_bins_vals),.(pos),summarise,meanrelcvg=mean(....relativeCvg),samples=paste0(unique(....sample),collapse=","))
    }
    if(format=="TCGA")
    {
     current_merged_df_bins_vals <- current_merged_df[,c("pos","current_gr.....Segment_Mean","....uuid")] %>% na.omit() %>% dplyr::mutate(current_gr.....Segment_Mean = current_gr.....Segment_Mean %>% as.character() %>% as.numeric())
     current_merged_df_bins_aggregated<-plyr::ddply(na.omit(current_merged_df_bins_vals),.(pos),summarise,segmentmean=mean(current_gr.....Segment_Mean),samples=paste0(unique(....uuid),collapse=","))
    }
    #current_merged_df_bins_aggregated_test<-ddply(na.omit(current_merged_df_bins_vals[1,]),.(pos),summarise,meanrelcvg=mean(current_merged_df_bins_vals$....relativeCvg))#,samples=list(unique(current_merged_df_bins_vals$....sample))
    
    #
    #insert bins that are not represented.
    unused_bins<-bins_underscored[!(bins_underscored %in% current_merged_df_bins_aggregated$pos)]
    unused_bins_rows<-as.data.frame(cbind(unused_bins,rep(0,length(unused_bins)),rep(samples[s],length(unused_bins))))
    if(format=="TARGET"){colnames(unused_bins_rows)<-c("pos","meanrelcvg","samples")
    unused_bins_rows$meanrelcvg<-as.numeric(unused_bins_rows$meanrelcvg)
    current_merged_df_bins_aggregated_with_unused<-rbind(current_merged_df_bins_aggregated[,c("pos","meanrelcvg","samples")],unused_bins_rows[,c("pos","meanrelcvg","samples")])
    
    }
    if(format=="TCGA"){colnames(unused_bins_rows)<-c("pos","segmentmean","samples")
    unused_bins_rows$segmentmean<-as.numeric(unused_bins_rows$segmentmean)
    current_merged_df_bins_aggregated_with_unused<-rbind(current_merged_df_bins_aggregated[,c("pos","segmentmean","samples")],unused_bins_rows[,c("pos","segmentmean","samples")])
    }
    current_merged_df_bins_aggregated_with_unused<-current_merged_df_bins_aggregated_with_unused[order(underscored_pos_to_GRanges(current_merged_df_bins_aggregated_with_unused$pos)),]
    #end testing
    if(format=="TARGET"){
      
    relcvg_df<-as.data.frame(current_merged_df_bins_aggregated_with_unused$meanrelcvg)
    rownames(relcvg_df)<-current_merged_df_bins_aggregated_with_unused$pos
    colnames(relcvg_df)<-samples[s]
    print(paste0(samples[s]," complete"))
    return(relcvg_df)
    }
    if(format=="TCGA"){
      
      segmentmean_df<-as.data.frame(current_merged_df_bins_aggregated_with_unused$segmentmean)
      rownames(segmentmean_df)<-current_merged_df_bins_aggregated_with_unused$pos
      colnames(segmentmean_df)<-samples[s]
      print(paste0(samples[s]," complete"))
      return(segmentmean_df)
    }
    
    

  }
  return(sample_aggregated_segvals)
}
