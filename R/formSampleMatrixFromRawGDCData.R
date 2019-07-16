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
#' @importFrom stringr str_detect
#' @importFrom plyr ddply
#' @param tcga_files GDC files to be read
#' @param format file format, TCGA or TARGET.
#' @param binsize the binsize, in base pairs (default 1Mb or 1e6).  This value provides a good balance of resolution and speed with memory sensitive applications.
#' @param freadskip the number of lines to skip in the GDC files, typically 14 (the first 13 lines are metadata and the first is a blank line in NBL data). Adjust as needed.
#' @param debug debug mode enable (allows specific breakpoints to be checked).
#' @param chromosomes A vector of chromosomes to be used. Defaults to chr1-chrX,
#'  but others can be added e.g. chrY or chrM for Y chromosome or mitochondrial DNA.
#'   Format expected is a character vector, e.g. c("chr1", "chr2", "chr3").
#' @return sample_aggregated_segvals A dataframe containing the aggregated segmentation values, based on the parameters provided.
#' @examples 
#' #Pipeline examples would be too large to include in package checks.
#' #please see browseVignettes("CNVScope") for a demonstration.
#' 
#' @export
globalVariables(c('begin','s',".","pos",'....relativeCvg','....sample','current_gr.....Segment_Mean','....uuid'),add=F)

formSampleMatrixFromRawGDCData<-function(tcga_files=NULL,format="TARGET",binsize=1e6,freadskip=NULL, parallel = F,debug=F,chromosomes=paste0("chr",c(seq(1:22),"X"),"_"))
{
  #chromosomes<-paste0("chr",c(seq(1:22),"X"),"_")
  # TCGA_CNV_data_with_sample_info<-ldply(tcga_files,
  #                                       function(x) {input_csv<-fread(x,skip=freadskip)
  #                                       sample_info_colsplit<-reshape2::colsplit(basename(x),"_|-|\\.",c("pre","project","num","sample","comparison","fn_ext"))
  #                                       input_csv_with_sample_info<-dplyr::bind_cols(input_csv,sample_info_colsplit[rep(1,nrow(input_csv)),])
  #                                       return(input_csv_with_sample_info)
  #                                       }
  # )
  #creates chromosomes object. This is necessary to clean out 
  #non-biological chromosomes of the form   chrUn_gl000211

  if(format=="TARGET")
  {
    if(is.null(freadskip)) {freadskip=14}
    TCGA_CNV_data_with_sample_info<-plyr::ldply(tcga_files,function(x) freadGDCfile(x,fread_skip=freadskip)) 
    #TCGA_CNV_data_with_sample_info <- tcga_files %>% purrr::map_dfr(freadGDCfile, format = "TCGA",fread_skip=freadskip)
    #combines the files in a dataframe with ldply and a modified fread.
  }
  
  if(format=="TCGA"){
    TCGA_CNV_data_with_sample_info<-plyr::ldply(tcga_files,freadGDCfile,format="TCGA") #TCGA format files are really simple to read, compared to TARGET files.
    #TCGA_CNV_data_with_sample_info$....sample<-tcga_files
    #TCGA_CNV_data_with_sample_info <- tcga_files %>% purrr::map_dfr(freadGDCfile, format = "TCGA",fread_skip=freadskip)
    #return(TCGA_CNV_data_with_sample_info)
    #combines the files in a dataframe with ldply and a modified fread.
  }
  
  #TCGA_CNV_data_with_sample_info_small<-ldply(tcga_files[1:100],freadGDCfile)  
  #end testing
  TCGA_CNV_data<-TCGA_CNV_data_with_sample_info
  if(format=="TCGA"){
    if(!(all(stringr::str_detect(TCGA_CNV_data$Chromosome,"chr")))){
      TCGA_CNV_data$Chromosome<-paste0("chr",TCGA_CNV_data$Chromosome) #1->chr1, 9->chr9.
    }
    colnames(TCGA_CNV_data)<-gsub("Chromosome",">chr",gsub("End","end",gsub("Start","begin",colnames(TCGA_CNV_data)))) #adds >chr, begin, and end columns if they come in a different form.
  }
  if(debug){browser()}
  TCGA_CNV_data_range_filtered<-TCGA_CNV_data %>% tidyr::drop_na(begin,end) 
  #drops those with a missing begin and end.
  TCGA_CNV_data_dt<-data.table::as.data.table(TCGA_CNV_data_range_filtered)
  #converts to data table.
  TCGA_CNV_data_gr<-GenomicRanges::GRanges(seqnames = TCGA_CNV_data_range_filtered$`>chr`,ranges = IRanges::IRanges(start = TCGA_CNV_data_range_filtered$begin,end = TCGA_CNV_data_range_filtered$end),... = TCGA_CNV_data_range_filtered[,4:ncol(TCGA_CNV_data_range_filtered)])
  #creates GRanges object with other columns appended. These can be accessed using mcols()
  bins<-GenomicRanges::tileGenome(seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens),tilewidth=binsize,cut.last.tile.in.chrom = T)
  #creates bins using the tileGenome function.
  if(debug){browser()}
  bins<-bins[as.character(bins@seqnames) %in% gsub("_","",chromosomes)]
  #removes Y and junk chromosomes. #modify chromosomes object at the top to add a Y if needed.
  #Y is best removed unless ALL of the participants are male, e.g. prostate cancer.
  rownames_gr = bins
  colnames_gr = bins
  #sets row and column bins
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
  bins_underscored<-GRanges_to_underscored_pos(bins) #convert bins to underscored strings.
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
      #grabs the GRanges object for the sth uuid.
    } 
    if(format == "TARGET")
    {
      current_gr<-TCGA_CNV_data_gr[mcols(TCGA_CNV_data_gr)$....sample %in% samples[s]]
      #grabs the current GRanges object for TARGET data.
    }
    current_merged_df<-as.data.frame(IRanges::mergeByOverlaps(bins,current_gr))
    #takes the range and merges it over the bins.
    current_merged_df$pos<-unlist(
      tidyr::unite(
      current_merged_df[,c("bins.seqnames","bins.start","bins.end")],pos
      )
      ) #takes three columns of chr, start, end and converts it to "chr_1000_2000" format.
    #sort(table(current_merged_df$pos),decreasing=T)
    if(format=="TARGET")
    {
    current_merged_df_bins_vals<-current_merged_df[,c("pos","....relativeCvg","....sample")] #,"....comparison"
    #grabs the BIN position, copy number, and sample of the current merged DF.
    current_merged_df_bins_vals$....relativeCvg<-as.numeric(as.character(current_merged_df_bins_vals$....relativeCvg)) #converts relative coverage to a number.
    current_merged_df_bins_vals<-na.omit(current_merged_df_bins_vals) #removes NAs
    current_merged_df_bins_aggregated<-plyr::ddply(na.omit(current_merged_df_bins_vals),.(pos),summarise,meanrelcvg=mean(....relativeCvg),samples=paste0(unique(....sample),collapse=",")) 
    #removes NAs from the merged by overlaps df, grabs the position, mean copy number, and concatenated samples for each bin like "NAKKEF,PALHRL,ABCABC". #note that it removes nonunique samples from this new string  of samples for each bin.
    }
    if(format=="TCGA")
    {
     current_merged_df_bins_vals <- current_merged_df[,c("pos","current_gr.....Segment_Mean","....uuid")] %>% na.omit() %>% dplyr::mutate(current_gr.....Segment_Mean = current_gr.....Segment_Mean %>% as.character() %>% as.numeric())
     #grabs position, CN number, sample ID
     #drops NAs
     #converts the CN from factor to character to numeric.
     current_merged_df_bins_aggregated<-plyr::ddply(na.omit(current_merged_df_bins_vals),.(pos),summarise,segmentmean=mean(current_gr.....Segment_Mean),samples=paste0(unique(....uuid),collapse=","))
     #removes NAs from the merged by overlaps df, grabs the position, mean copy number, and concatenated samples for each bin like "NAKKEF,PALHRL,ABCABC". #note that it removes nonunique samples from this new string of samples for each bin.
    }
    #current_merged_df_bins_aggregated_test<-ddply(na.omit(current_merged_df_bins_vals[1,]),.(pos),summarise,meanrelcvg=mean(current_merged_df_bins_vals$....relativeCvg))#,samples=list(unique(current_merged_df_bins_vals$....sample))
    
    #
    #insert bins that are not represented.
    unused_bins<-bins_underscored[!(bins_underscored %in% current_merged_df_bins_aggregated$pos)]
    unused_bins_rows<-as.data.frame(cbind(unused_bins,rep(0,length(unused_bins)),rep(samples[s],length(unused_bins)))) #sets the unused bins to have zero CN.
    if(format=="TARGET"){colnames(unused_bins_rows)<-c("pos","meanrelcvg","samples")
    unused_bins_rows$meanrelcvg<-as.numeric(unused_bins_rows$meanrelcvg) 
    #converts CN to numeric.
    current_merged_df_bins_aggregated_with_unused<-rbind(current_merged_df_bins_aggregated[,c("pos","meanrelcvg","samples")],unused_bins_rows[,c("pos","meanrelcvg","samples")])
    #combines the empty bins with the rest.
    }
    if(format=="TCGA"){colnames(unused_bins_rows)<-c("pos","segmentmean","samples")
    unused_bins_rows$segmentmean<-as.numeric(unused_bins_rows$segmentmean)
    #converts CN to numeric.
    current_merged_df_bins_aggregated_with_unused<-rbind(current_merged_df_bins_aggregated[,c("pos","segmentmean","samples")],unused_bins_rows[,c("pos","segmentmean","samples")])
    #combines the empty bins with the rest.
    
    }
    current_merged_df_bins_aggregated_with_unused<-current_merged_df_bins_aggregated_with_unused[order(underscored_pos_to_GRanges(current_merged_df_bins_aggregated_with_unused$pos)),]
    #orders the combined bins once more.
    #end testing
    if(format=="TARGET"){
      
    relcvg_df<-as.data.frame(current_merged_df_bins_aggregated_with_unused$meanrelcvg)
    rownames(relcvg_df)<-current_merged_df_bins_aggregated_with_unused$pos
    colnames(relcvg_df)<-samples[s]
    print(paste0(samples[s]," complete")) 
    #prints completion.
    #takes this merged dataframe, for the given sample and returns it to the loop. 
    return(relcvg_df)
    }
    if(format=="TCGA"){
  #TCGA follows precisely the same procedure as target.    
      segmentmean_df<-as.data.frame(current_merged_df_bins_aggregated_with_unused$segmentmean)
      rownames(segmentmean_df)<-current_merged_df_bins_aggregated_with_unused$pos
      colnames(segmentmean_df)<-samples[s]
      print(paste0(samples[s]," complete"))
      return(segmentmean_df)
    }
    
    

  }
  return(sample_aggregated_segvals)
  #each round of the loop gives each CN mean for each bin in a single sample
  #foreach combines the samples into a dataframe with cbind.
  #the column names are the samples and the rownames are the bin positions.
}
