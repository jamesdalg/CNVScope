#' Server component of the CNVScope plotly shiny application.
#'
#'  Server function of the CNVScope shiny application. run with runCNVScopeShiny
#' @name CNVScopeserver 
#' @keywords CNV heatmap shiny plotly
#' @import shinycssloaders shinythemes visNetwork ggplot2 reshape2 magrittr htmltools htmlwidgets jointseg logging foreach GenomicInteractions shinythemes
#' @importFrom tidyr unite
#' @rawNamespace import(circlize, except = degree)
#' @rawNamespace import(shiny, except = runExample)
#' @rawNamespace import(shinyjs, except = runExample)
#' @rawNamespace import(RCurl, except = reset)
#' @rawNamespace import(plotly, except = c(last_plot,select,filter))
#' @rawNamespace import(igraph, except = c(decompose, spectrum, groups))
#' @rawNamespace import(data.table, except = c(melt, dcast))
#' @rawNamespace import(GenomicFeatures ,except = show)
#' @param session The shiny session object for the application.
#' @param input shiny server input
#' @param output shiny server output
#' @param debug enable debugging mode
#' @return None

#' @examples
#' \dontrun{
#' runCNVScopeShiny()
#' }
#' @export
#globalVariables(c("ensembl_gene_tx_data_gr","baseurl","chromosomes","downsample_factor","basefn",
#                  "subset_name",
#                  "expression_data_gr_nbl",'start2','start1','value','Var1','Var2','value1',
#                 'tcga_type','census_data_gr','common_coords','myReactives',
#                  'genev','delete.isolates','freq_data'),add = F)
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."), add=F)
CNVScopeserver<-function(session,input, output, debug=F) {
ensembl_gene_tx_data_gr <- if(exists("ensembl_gene_tx_data_gr")){get("ensembl_gene_tx_data_gr")} else {NULL}
baseurl <- if(exists("baseurl")){get("baseurl")} else {NULL}
basefn <- if(exists("basefn")){get("basefn")} else {NULL}
osteofn <- if(exists("osteofn")){get("osteofn")} else {NULL}
start1 <- if(exists("start1")){get("start1")} else {NULL}
start2 <- if(exists("start2")){get("start2")} else {NULL}
value <- if(exists("value")){get("value")} else {NULL}
value1 <- if(exists("value1")){get("value1")} else {NULL}
Var1 <- if(exists("Var1")){get("Var1")} else {NULL}
Var2 <- if(exists("Var2")){get("Var2")} else {NULL}
bins.seqnames <- if(exists("bins.seqnames")){get("bins.seqnames")} else {NULL}
bins.start <- if(exists("bins.start")){get("bins.start")} else {NULL}
bins.end <- if(exists("bins.end")){get("bins.end")} else {NULL}
expression_data_gr <- if(exists("expression_data_gr")){get("expression_data_gr")} else {NULL}
common_coords <- if(exists("common_coords")){get("common_coords")} else {NULL}
myReactives <- if(exists("myReactives")){get("myReactives")} else {NULL}
genev <- if(exists("genev")){get("genev")} else {NULL}
delete.isolates <- function(graph, mode = 'all') {
  isolates <- which(igraph::degree(graph, mode = mode) == 0) 
  delete.vertices(graph, isolates)
}
freq_data <- if(exists("freq_data")){get("freq_data")} else {NULL}
#adjpvalue chr cn correlation genes_text probe visval
adjpvalue <- if(exists("adjpvalue")){get("adjpvalue")} else {NULL}
chr <- if(exists("chr")){get("chr")} else {NULL}
cn <- if(exists("cn")){get("cn")} else {NULL}
correlation <- if(exists("correlation")){get("correlation")} else {NULL}
genes_text <- if(exists("genes_text")){get("genes_text")} else {NULL}
probe <- if(exists("probe")){get("probe")} else {NULL}
visval <- if(exists("visval")){get("visval")} else {NULL}
  privpolurl <- a("NCI Privacy Policy", href="https://www.cancer.gov/policies/privacy-security",target="_blank")
  output$privpol <- renderUI({
    tagList(privpolurl)})
  downsample_factor<-NULL
  subset_name<-NULL
  #expression_data_gr_nbl<-NULL
  tcga_type<-NULL
  chrom.pairs<-NULL
  
  
  printLogJs <- function(x, ...) {
    
    logjs(x)
    
    T
  }
  observe({
    if (input$geneSearch == 0) {return()}
    
    x<-isolate(input$geneSearch)
    
    #browser()
    if(x!=0 & isolate(input$gene_input_col)!=""& isolate(input$gene_input_row)!=""){
      if(length(ensembl_gene_tx_data_gr[ensembl_gene_tx_data_gr$....external_gene_name==isolate(input$gene_input_col)])!=0) {
        
        colgene_loc<-paste0(ensembl_gene_tx_data_gr[ensembl_gene_tx_data_gr$....external_gene_name==isolate(input$gene_input_col)][1]$....chromosome_name,":",
                            ensembl_gene_tx_data_gr[ensembl_gene_tx_data_gr$....external_gene_name==isolate(input$gene_input_col)][1]$....start_position,"-",
                            ensembl_gene_tx_data_gr[ensembl_gene_tx_data_gr$....external_gene_name==isolate(input$gene_input_col)][1]$....end_position)
      } else {
        
        colgene_loc<-""}
      
      if(length(ensembl_gene_tx_data_gr[ensembl_gene_tx_data_gr$....external_gene_name==isolate(input$gene_input_row)])!=0) {
        
        
        rowgene_loc<-paste0(ensembl_gene_tx_data_gr[ensembl_gene_tx_data_gr$....external_gene_name==isolate(input$gene_input_row)][1]$....chromosome_name,":",
                            ensembl_gene_tx_data_gr[ensembl_gene_tx_data_gr$....external_gene_name==isolate(input$gene_input_row)][1]$....start_position,"-",
                            ensembl_gene_tx_data_gr[ensembl_gene_tx_data_gr$....external_gene_name==isolate(input$gene_input_row)][1]$....end_position)
      } else {
        #   #browser()
        rowgene_loc<-""}
      updateTextInput(session,"loc_input_col",value=colgene_loc)
      updateTextInput(session,"loc_input_row",value=rowgene_loc)
      if(length(ensembl_gene_tx_data_gr[ensembl_gene_tx_data_gr$....external_gene_name==isolate(input$gene_input_row)])!=0){
        updateSelectInput(session,"chrom2",selected = paste0(ensembl_gene_tx_data_gr[ensembl_gene_tx_data_gr$....external_gene_name==isolate(input$gene_input_row)][1]$....chromosome_name,"_"))}
      if(length(ensembl_gene_tx_data_gr[ensembl_gene_tx_data_gr$....external_gene_name==isolate(input$gene_input_col)])!=0){
        updateSelectInput(session,"chrom1",selected = paste0(ensembl_gene_tx_data_gr[ensembl_gene_tx_data_gr$....external_gene_name==isolate(input$gene_input_col)][1]$....chromosome_name,"_"))}
      
    } #end check to see if there is input in the gene search.
    
  })
  observeEvent(event_data("plotly_click"), {
    showTab(inputId = "tabs",select = T, target = "sample info")
    if(isolate(input$data_source)=="linreg_osteosarcoma_CNVkit")
    {
      
      showTab(inputId="tabs",target="gain/loss frequency")
    }
    showTab(inputId="tabs",target="sample info")
    showTab(inputId="tabs",target="COSMIC cancer gene census") 
    showTab(inputId="tabs",target="expression_data")
  })
  observeEvent(input$goButton, {
    showTab(inputId = "tabs",select = T, target = "Plots")
    if(isolate(input$data_source)!="linreg_osteosarcoma_CNVkit")
    {
      hideTab(inputId="tabs",target="gain/loss frequency")
    }
    if(isolate(input$data_source)=="linreg_osteosarcoma_CNVkit")
    {
      showTab(inputId="tabs",select = F,target="gain/loss frequency")
    }
  })
  observeEvent(input$data_source, {
    if(isolate(input$data_source)!="linreg_osteosarcoma_CNVkit")
    {
      hideTab(inputId="tabs",target="gain/loss frequency")
    }
    if(isolate(input$data_source)=="linreg_osteosarcoma_CNVkit")
    {
      showTab(inputId="tabs",select = F,target="gain/loss frequency")
    }
    if(is.null(event_data("plotly_click"))){
      hideTab(inputId="tabs",target="gain/loss frequency")
      hideTab(inputId="tabs",target="sample info")
      hideTab(inputId="tabs",target="COSMIC cancer gene census") 
      hideTab(inputId="tabs",target="expression_data")
    }
  })
  getHeight<-function()
  {
    return(isolate(input$heatmapHeight)) 
  }
  addHandler(printLogJs)
  isolate(input$goButton)
  # observe({
  #   input$goButton
  #   if(!is.null(isolate(input$loc_input_row))){
  # updateSelectInput(session,"chrom1",chromosomes,selected=paste0(as.character(GRanges(isolate(input$loc_input_row))@seqnames),"_"))}
  # })
  output$plotlyChromosomalHeatmap <- renderPlotly({
    
    if (input$goButton == 0) {return()}
    
    input$goButton
    # if(!file.exists(
    #   (
    #     paste0(getwd(),"/matrix/linreg/",
    #            chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom1))))],chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom2))))],
    #            "melted_downsampled_linreg.RData")
    #   )
    #   )){ return("file does not exist!");}
    #if there is location data, change the chromosomes from what they were chosen.
    # 
    #isolate(input$loc_input_row)
    # observe({
    # updateSelectInput(session,"chrom1",chromosomes,selected=paste0(as.character(GRanges(isolate(input$loc_input_row))@seqnames),"_"))
    # })
    
    
    if(isolate(input$data_source)=="linreg_osteosarcoma_CNVkit")
    {     load( url(paste0(paste0(baseurl,"matrix/linreg/unrescaled/",
                                  chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom1))))],chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom2))))],
                                  "melted_downsampled_linreg_unrescaled.RData"))))
      
      load( url(paste0(paste0(baseurl,"matrix/linreg/unrescaled/full/",
                              chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom1))))],chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom2))))],
                              "melted_full_linreg_max_cap_75.RData"))))
      #browser()
      downsample_factor<<-4
      tryCatch(bin_data<-readRDS((url(paste0(baseurl,"bin_data.rds")))),error = function(e) NULL) 
      tryCatch(bin_data<-readRDS((paste0(osteofn,"bin_data.rds"))),error = function(e) NULL) 
      
    }
    # 
    if(isolate(input$data_source)=="TCGA_SARC_SNP6")
    {
      load( url(paste0(paste0(baseurl,"matrix/TCGA_SARC/downsampled_factor_8/",
                              chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom1))))],chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom2))))],
                              "melted_downsampled_TGCA_SARC_unrescaledv2.RData"))))
      
      # load( url(paste0(paste0(baseurl,"matrix/TCGA_SARC/full/",
      #                  chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom1))))],chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom2))))],
      #                  "melted_full_TGCA_SARC_unrescaled.RData"))))
      downsample_factor<<-8
    }
    #
    if(isolate(input$data_source)=="TCGA_BRCA_low_pass")
    {
      #
      sample_name<-"BRCA_output_matrix1e6"
      load( url(paste0(paste0(baseurl,"matrix/TCGA_low_pass/BRCA/",
                              paste0(chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom1))))],chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom2))))],"melted_downsampled_TGCA_",sample_name,"_unrescaled",".RData")
      ))))
      ggplotmatrix_full<-ggplotmatrix
    }
    if(isolate(input$data_source)=="TCGA_AML_low_pass")
    {
      sample_name<-"AML_output_matrix1e6"
      load( url(paste0(paste0(baseurl,"matrix/TCGA_low_pass/AML/",
                              paste0(chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom1))))],chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom2))))],"melted_downsampled_TGCA_",sample_name,"_unrescaled",".RData")
      ))))
      ggplotmatrix_full<-ggplotmatrix
    }
    if(isolate(input$data_source)=="TCGA_PRAD_low_pass")
    {
      sample_name<-"PRAD_output_matrix1e6"
      load( url(paste0(paste0(baseurl,"matrix/TCGA_low_pass/PRAD/",
                              paste0(chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom1))))],chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom2))))],"melted_downsampled_TGCA_",sample_name,"_unrescaled",".RData")
      ))))
      ggplotmatrix_full<-ggplotmatrix
    }
    if(isolate(input$data_source)=="TCGA_NBL_low_pass")
    {
      sample_name<-"NBL_output_matrix1e6"
      load( paste0(paste0(basefn,"matrix/TCGA_low_pass/NBL/",
                              paste0(isolate(input$chrom1),isolate(input$chrom2),"nbl_sample_matched_unrescaled.RData")
      )))
      #browser()
      #     ggplotmatrix
      ggplotmatrix_full<-ggplotmatrix
      tryCatch(bin_data<<-readRDS((url(paste0(baseurl,"bin_data_nbl.rds")))),error = function(e) NULL) 
      tryCatch(bin_data<<-readRDS((paste0(basefn,"bin_data_nbl.rds"))),error = function(e) NULL) 
      
    }
    if(isolate(input$data_source) %in% c("TCGA_NBL_stage3_subset","TCGA_NBL_stage4_subset","TCGA_NBL_stage4s_subset","TCGA_NBL_myc_amp_subset","TCGA_NBL_not_myc_amp_subset"))
    {
      #browser()
      subset_name<<-gsub("_subset","",gsub("TCGA_NBL_","",paste0(input$data_source)))
      sample_name<-"NBL_output_matrix1e6"
      # load( url(paste0(paste0(baseurl,"matrix/TCGA_low_pass/NBL/",subset_name,"/",
      #                         paste0(chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom1))))],chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom2))))],"melted_downsampled_TGCA_","NBLsample_matched","_unrescaled",subset_name,".RData")
      # ))))
      load( url(paste0(paste0(baseurl,"matrix/TCGA_low_pass/NBL/",subset_name,"/",
                              paste0(isolate(input$chrom1),isolate(input$chrom2),"melted_downsampled_TGCA_","NBLsample_matched","_unrescaled",subset_name,"pos_neg.RData")
      ))))
      if(length(bin_data$probe)==0)
      {
        bin_data$probe<-rownames(bin_data)
      }
      
      
      
      
      
      #      ggplotmatrix
      ggplotmatrix_full<-ggplotmatrix
      tryCatch(bin_data<<-readRDS((url(paste0(baseurl,"bin_data_nbl_",subset_name,".rds")))),error = function(e) NULL) 
      tryCatch(bin_data<<-readRDS((paste0(basefn,"bin_data_nbl_",subset_name,".rds"))),error = function(e) NULL) 
      input_mat<-bin_data %>% dplyr::select(-probe)
      rownames(input_mat)<-bin_data$probe
      #
      tryCatch(expression_data_gr_nbl<<-readRDS(url(paste0(baseurl,"tcga_nbl_expression_",subset_name,"subset.rds"))),error = function(e) NULL)
      tryCatch(expression_data_gr_nbl<<-readRDS(paste0(basefn,"tcga_nbl_expression_",subset_name,"subset.rds")),error = function(e) NULL)
      
      #server-side processing(disabled):
      # tryCatch(tcga_gr<<-readRDS((url(paste0(baseurl,"tcga_gr_no_stats.rds")))),error = function(e) NULL) 
      # tryCatch(tcga_gr<<-readRDS((paste0(basefn,"tcga_gr_no_stats.rds"))),error = function(e) NULL) 
      # tryCatch(tcga_dfs_cbind_with_ensg_with_ensembl_fpkm<<-readRDS((url(paste0(baseurl,"tcga_dfs_cbind_with_ensg_with_ensembl_fpkm_caseid.rds")))),error = function(e) NULL) 
      # tryCatch(tcga_dfs_cbind_with_ensg_with_ensembl_fpkm<<-readRDS((paste0(basefn,"tcga_dfs_cbind_with_ensg_with_ensembl_fpkm_caseid.rds"))),error = function(e) NULL) 
      # 
      # tcga_dfs_cbind_with_ensg_with_ensembl_fpkm_subset<-as.data.frame(tcga_dfs_cbind_with_ensg_with_ensembl_fpkm)[,na.omit(match(colnames(bin_data),colnames(tcga_dfs_cbind_with_ensg_with_ensembl_fpkm)))]
      # #dim(tcga_dfs_cbind_with_ensg_with_ensembl_fpkm_subset)
      # mcols(tcga_gr)$rowMean<-rowMeans(tcga_dfs_cbind_with_ensg_with_ensembl_fpkm_subset) #tcga_dfs_cbind_with_ensg[,2:ncol(tcga_dfs_cbind_with_ensg)]
      # mcols(tcga_gr)$rowMeanPctl<-heatmaply::percentize(rowMeans(tcga_dfs_cbind_with_ensg_with_ensembl_fpkm_subset))
      # mcols(tcga_gr)$rowVar<-matrixStats::rowVars(as.matrix(tcga_dfs_cbind_with_ensg_with_ensembl_fpkm_subset))
      # mcols(tcga_gr)$rowVarPctl<-heatmaply::percentize(matrixStats::rowVars(as.matrix(tcga_dfs_cbind_with_ensg_with_ensembl_fpkm_subset)))
      # mcols(tcga_gr)$SYMBOL<-mcols(tcga_gr)$....external_gene_name
      # mcols(tcga_gr)$gene_type<-mcols(tcga_gr)$....gene_biotype
      # expression_data_gr<<-tcga_gr
    }
    if(isolate(input$data_source)=="TCGA_OS_low_pass")
    {
      sample_name<-"OS_output_matrix1e6"
      load( url(paste0(paste0(baseurl,"matrix/TCGA_low_pass/OS/",
                              paste0(chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom1))))],chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom2))))],"melted_downsampled_TGCA_",sample_name,"_unrescaled",".RData")
      ))))
      ggplotmatrix_full<-ggplotmatrix
    }
    
    ggplotmatrix$value<-signedRescale(ggplotmatrix$value,max_cap=isolate(input$max_cap))[,1]
    ggplotmatrix<-dplyr::bind_cols(ggplotmatrix,reshape2::colsplit(ggplotmatrix$Var1,"_",c("chr1","start1","end1")))
    ggplotmatrix<-dplyr::bind_cols(ggplotmatrix,reshape2::colsplit(ggplotmatrix$Var2,"_",c("chr2","start2","end2")))
    ggplotmatrix<-ggplotmatrix[order(ggplotmatrix$start1,ggplotmatrix$start2),]
    if(!is.null(ggplotmatrix)){ggplotmatrix<<-ggplotmatrix}
    #
    if(!is.null(ggplotmatrix_full)){ ggplotmatrix_full$value<-signedRescale(ggplotmatrix_full$value,max_cap=isolate(input$max_cap))[,1]}
    if(!is.null(ggplotmatrix_full)){ggplotmatrix_full<<-ggplotmatrix_full}
    recast_matrix<-reshape2::dcast(data=ggplotmatrix,formula=Var1 ~ Var2, var = ggplotmatrix$value) #this creates a matrix in wide format.
    if(ncol(recast_matrix)!=nrow(recast_matrix))
    {
      rownames(recast_matrix)<-recast_matrix$Var1
      recast_matrix<-recast_matrix[,2:ncol(recast_matrix)]
    }
    #
    recast_matrix_full<-reshape2::dcast(data=ggplotmatrix_full,formula=Var1 ~ Var2, var = ggplotmatrix_full$value) #this creates a matrix with 
    if(ncol(recast_matrix_full)!=nrow(recast_matrix_full))
    {
      rownames(recast_matrix_full)<-recast_matrix_full$Var1
      recast_matrix_full<-recast_matrix_full[,2:ncol(recast_matrix_full)]
    }
    #
    #resorting recast_matrix
    if(!is.null(recast_matrix)){recast_matrix<<-recast_matrix}
    if(!is.null(recast_matrix_full)){recast_matrix_full<<-recast_matrix_full}
    rownames_gr<-underscored_pos_to_GRanges(rownames(recast_matrix),zeroToOneBasedStart = F,zeroToOneBasedEnd = F)
    colnames_gr<-underscored_pos_to_GRanges(colnames(recast_matrix),zeroToOneBasedStart = F,zeroToOneBasedEnd = F)
    rownames_gr_full<-underscored_pos_to_GRanges(rownames(recast_matrix_full),zeroToOneBasedStart = F,zeroToOneBasedEnd = F)
    colnames_gr_full<-underscored_pos_to_GRanges(colnames(recast_matrix_full),zeroToOneBasedStart = F,zeroToOneBasedEnd = F)
    if(!is.null(rownames_gr)){rownames_gr<<-rownames_gr}
    if(!is.null(rownames_gr_full)){rownames_gr_full<<-rownames_gr_full}
    if(!is.null(colnames_gr)){colnames_gr<<-colnames_gr}
    if(!is.null(colnames_gr_full)){colnames_gr_full<<-colnames_gr_full}
    
    ggplotmatrix$value1<-gsub("col genes:","row genes:",ggplotmatrix$value1)
    ggplotmatrix$value1<-gsub("row_genes:","col_genes:",ggplotmatrix$value1)
    rownames_ordered<-GRanges_to_underscored_pos(rownames_gr[order(rownames_gr)])
    colnames_ordered<-GRanges_to_underscored_pos(colnames_gr[order(colnames_gr)])
   if(debug){browser()}
    recast_matrix<-recast_matrix[rownames_ordered,colnames_ordered]
    block_indices_row<-jointseg::jointSeg(recast_matrix,K=10,method="RBS")$bestBkp
    block_indices_col<-jointseg::jointSeg(t(recast_matrix),K=10,method="RBS")$bestBkp
    block_index_labels_row<-rownames(recast_matrix)[block_indices_row]
    block_index_labels_col<-colnames(recast_matrix)[block_indices_col]
    
    # xfactor<-as.factor(ggplotmatrix$Var1)
    # levels(xfactor)<-order(colnames_gr)
    # yfactor<-as.factor(ggplotmatrix$Var1)
    # levels(yfactor)<-order(rownames_gr)
    # p <- ggplot(data = ggplotmatrix ) + #geom_tile() + theme_void()
    #   geom_raster(aes(x = xfactor, y = yfactor,fill=value,text=paste0("value:",value,"\nrow:",Var1,"\ncol:",Var2,"\n",value1))) + scale_x_discrete(breaks = block_index_labels_col) +
    #   scale_y_discrete(breaks = block_index_labels_row) + theme(axis.text.x = element_text(angle=60, hjust=1)) +  
    #   ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.5, limits = c(0, 1)) +  theme(legend.position="bottom",axis.title = element_blank()) + coord_flip() #+ scale_y_reverse(breaks=block_indices)
    # 

#recreate input matrix, add rownames.
    options(stringsAsFactors = F)
input_mat<-bin_data %>% dplyr::select(-probe) %>% as.data.frame()
rownames(input_mat)<-bin_data$probe
#correlate input matrix
if(isolate(input$cor_method)!="spearman - pearson"){
input_mat_cor<-cor(t(input_mat),method=isolate(input$cor_method))
} else {
  input_mat_cor<-cor(t(input_mat),method="spearman")-cor(t(input_mat),method="pearson")
}
#wide to long
input_mat_cor_flat<-input_mat_cor %>% reshape2::melt()
#grab ggplotmatrix and add correlation values.
#if(!isolate(input$genes_toggle)){ggplotmatrix$value1<-NULL}
ggplotmatrix_joined<- dplyr::inner_join(x=ggplotmatrix,y=input_mat_cor_flat,by=c("Var1"="Var1","Var2"="Var2"))
colnames(ggplotmatrix_joined) <- ggplotmatrix_joined %>% colnames() %>%
  gsub(pattern = "value.x",replacement = "linregval") %>%
  gsub(pattern = "value.y",replacement = "correlation")
#convert the negative log p-values to p-values and apply two kinds of FDR correction.

ggplotmatrix_joined$pvalue<-exp(-(abs(ggplotmatrix_joined$orig_value)))
ggplotmatrix_joined$adjpvaluechr<-p.adjust(p = ggplotmatrix_joined$pvalue,method = "fdr")
ggplotmatrix_joined$adjpvaluegenome<-p.adjust(p = ggplotmatrix_joined$pvalue,method = "fdr",
                                              n = dim(input_mat)[1]*dim(input_mat)[2])
ggplotmatrix_joined<<-ggplotmatrix_joined
rownames_ordered<-GRanges_to_underscored_pos(rownames_gr[order(rownames_gr)])
colnames_ordered<-GRanges_to_underscored_pos(colnames_gr[order(colnames_gr)])
if(isolate(input$fdr_correction)=="chromosome_pair"){
  ggplotmatrix_joined$adjpvalue<-ggplotmatrix_joined$adjpvaluechr    
} else {
  if(isolate(input$fdr_correction)=="genome"){
ggplotmatrix_joined$adjpvalue<-ggplotmatrix_joined$adjpvaluegenome  
}
}
ggplotmatrix_joined<<-ggplotmatrix_joined
if(isolate(input$visval)=="Correlation") {
  ggplotmatrix_joined$visval<-ggplotmatrix_joined$correlation
} else {
  if(isolate(input$visval)=="-log(Linear Regression P-value) * correlation sign") {
  ggplotmatrix_joined$visval<-ggplotmatrix_joined$linregval
  }
}
if(isolate(input$pval_filter_toggle)){
ggplotmatrix_joined$visval<-ifelse(ggplotmatrix_joined$adjpvalue<0.05,ggplotmatrix_joined$linregval,0.5)
} else {
  ggplotmatrix_joined$visval<-ggplotmatrix_joined$linregval
}
if(!isolate(input$genes_toggle)){
  ggplotmatrix_joined$genes_text<-rep("",nrow(ggplotmatrix_joined))
} else {
  ggplotmatrix_joined$genes_text<-ggplotmatrix_joined$value1
}
    #as.integer(as.character(reshape2::colsplit(ggplotmatrix$Var2,"_",c("chr2","start2","end2"))$start2))
    p <- ggplot(data = ggplotmatrix_joined ) + #geom_tile() + theme_void()
      geom_tile(aes(x =      as.numeric(start2),
                      y =      as.numeric(start1),
                      fill=visval,text=paste0("value:",visval,"\nrow:",Var1,"\ncol:",Var2,"\n",genes_text,"\nFDR p=",adjpvalue,"\n",isolate(input$cor_method)," Correlation=",correlation)),alpha=ifelse(ggplotmatrix_joined$adjpvaluechr<0.05,1.0,0.1)) + #
      scale_x_continuous(breaks = reshape2::colsplit(block_index_labels_col,"_",c("chr","start","end"))$start,labels = block_index_labels_col) +
      scale_y_continuous(breaks = reshape2::colsplit(block_index_labels_row,"_",c("chr","start","end"))$start,labels = block_index_labels_row) + theme(axis.text.x = element_text(angle=60, hjust=1)) +  
      ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.5, limits = c(0, 1)) +  theme(legend.position="bottom",axis.title = element_blank()) #+ geom_contour(binwidth = .395,aes(z=value))
###    browser()
    #+ coord_flip() #+ scale_y_reverse(breaks=block_indices)
    #p
    #lumpy_points_toggle
    if(isolate(input$data_source)=="linreg_osteosarcoma_CNVkit")
    {
      if(exists("osteofn"))
      {
        tryCatch(SVs_data_in_submatrix_coords<-readRDS(paste0(osteofn,"breakpoint_gint/",isolate(input$chrom1),isolate(input$chrom2),"SVs_data_in_submatrix_coords.rds" )),error = function(e) NULL) 
        tryCatch(lumpy_summarized_counts<-readRDS(paste0(osteofn,"lumpy_sv/",gsub("_","",isolate(input$chrom1)),gsub("_","",isolate(input$chrom2)),"SVs_data_in_submatrix_coords_lumpy_mirror.rds" )),error = function(e) NULL)    
      }else {
        tryCatch(SVs_data_in_submatrix_coords<-readRDS(url(paste0(baseurl,"breakpoint_gint/",isolate(input$chrom1),isolate(input$chrom2),"SVs_data_in_submatrix_coords.rds" ))),error = function(e) NULL) 
        tryCatch(lumpy_summarized_counts<-readRDS(url(paste0(baseurl,"lumpy_sv/",gsub("_","",isolate(input$chrom1)),gsub("_","",isolate(input$chrom2)),"SVs_data_in_submatrix_coords_lumpy_mirror.rds" ))),error = function(e) NULL)   
      }
      
    }
    if(isolate(input$data_source) %in% c("TCGA_AML_low_pass","TCGA_BRCA_low_pass","TCGA_OS_low_pass","TCGA_NBL_low_pass","TCGA_PRAD_low_pass"))
    {
      if(exists("basefn"))
      {
        tryCatch(SVs_data_in_submatrix_coords<-readRDS(paste0(basefn,"breakpoint_gint/TCGA_low_pass/",isolate(input$chrom1),isolate(input$chrom2),"SVs_data_in_submatrix_coords_common_coords.rds" )),error = function(e) NULL)
        tryCatch(lumpy_summarized_counts<-readRDS(paste0(basefn,"lumpy_sv/TCGA_low_pass/",isolate(input$chrom1),isolate(input$chrom2),"SVs_data_in_submatrix_coords_lumpy_mirror_TCGA_common_coords.rds" )),error = function(e) NULL)
        tcga_type<<-gsub("_low_pass","",gsub("TCGA_","",isolate(input$data_source)))
        tryCatch(TCGA_low_pass_sample_info<<-readRDS(paste0(basefn,"sample_info/",tcga_type,"TCGA_merged_dtv2.rds" )),error = function(e) NULL)
        if(exists("TCGA_low_pass_sample_info")){TCGA_low_pass_sample_info$pos<- tidyr::unite(TCGA_low_pass_sample_info,pos,bins.seqnames,bins.start,bins.end)$pos}
      } else {
        tryCatch(SVs_data_in_submatrix_coords<-readRDS(url(paste0(baseurl,"breakpoint_gint/TCGA_low_pass/",isolate(input$chrom1),isolate(input$chrom2),"SVs_data_in_submatrix_coords_common_coords.rds" ))),error = function(e) NULL)
        tryCatch(lumpy_summarized_counts<-readRDS(url(paste0(baseurl,"lumpy_sv/TCGA_low_pass/",isolate(input$chrom1),isolate(input$chrom2),"SVs_data_in_submatrix_coords_lumpy_mirror_TCGA_common_coords.rds" ))),error = function(e) NULL)
        tcga_type<<-gsub("_low_pass","",gsub("TCGA_","",isolate(input$data_source)))
        tryCatch(TCGA_low_pass_sample_info<<-readRDS(url(paste0(baseurl,"sample_info/",tcga_type,"TCGA_merged_dtv2.rds" ))),error = function(e) NULL)
        if(exists("TCGA_low_pass_sample_info")){TCGA_low_pass_sample_info$pos<- tidyr::unite(TCGA_low_pass_sample_info,pos,bins.seqnames,bins.start,bins.end)$pos}
      }
    }
    if(isolate(input$data_source) %in% c("TCGA_NBL_stage3_subset","TCGA_NBL_stage4_subset","TCGA_NBL_stage4s_subset","TCGA_NBL_myc_amp_subset","TCGA_NBL_not_myc_amp_subset"))
    {
      subset_name<<-gsub("_subset","",gsub("TCGA_NBL_","",paste0(input$data_source)))
      if(exists("basefn"))
      {
        tryCatch(SVs_data_in_submatrix_coords<-readRDS(paste0(basefn,"breakpoint_gint/TCGA_low_pass/",isolate(input$chrom1),isolate(input$chrom2),"SVs_data_in_submatrix_coords_common_coords.rds" )),error = function(e) NULL)
        tryCatch(lumpy_summarized_counts<-readRDS(paste0(basefn,"lumpy_sv/TCGA_low_pass/",isolate(input$chrom1),isolate(input$chrom2),"SVs_data_in_submatrix_coords_lumpy_mirror_TCGA_common_coords.rds" )),error = function(e) NULL)
        tcga_type<<-gsub("_low_pass","",gsub("TCGA_","",isolate(input$data_source)))
        tryCatch(TCGA_low_pass_sample_info<<-readRDS(paste0(basefn,"sample_info/",tcga_type,"TCGA_merged_dtv2.rds" )),error = function(e) NULL)
        if(exists("TCGA_low_pass_sample_info")){TCGA_low_pass_sample_info$pos<- tidyr::unite(TCGA_low_pass_sample_info,pos,bins.seqnames,bins.start,bins.end)$pos}
      } else {
        tryCatch(SVs_data_in_submatrix_coords<-readRDS(url(paste0(baseurl,"breakpoint_gint/TCGA_low_pass/",isolate(input$chrom1),isolate(input$chrom2),"SVs_data_in_submatrix_coords_common_coords.rds" ))),error = function(e) NULL)
        tryCatch(lumpy_summarized_counts<-readRDS(url(paste0(baseurl,"lumpy_sv/TCGA_low_pass/",isolate(input$chrom1),isolate(input$chrom2),"SVs_data_in_submatrix_coords_lumpy_mirror_TCGA_common_coords.rds" ))),error = function(e) NULL)
        tcga_type<<-gsub("_low_pass","",gsub("TCGA_","",isolate(input$data_source)))
        tryCatch(TCGA_low_pass_sample_info<<-readRDS(url(paste0(baseurl,"sample_info/",tcga_type,"TCGA_merged_dtv2.rds" ))),error = function(e) NULL)
        if(exists("TCGA_low_pass_sample_info")){TCGA_low_pass_sample_info$pos<- tidyr::unite(TCGA_low_pass_sample_info,pos,bins.seqnames,bins.start,bins.end)$pos}
      }
    }
    
    # return(lumpy_summarized_counts)
    #}
    
    #DISABLING CLIENT SIDE PROCESSING OF GenomicInteraction data.
    # submat_row_gr<-underscored_pos_to_GRanges(rownames(recast_matrix))
    # submat_col_gr<-underscored_pos_to_GRanges(colnames(recast_matrix))
    # breakpoint_gint_full_subset<-breakpoint_gint_full[anchorOne(breakpoint_gint_full)@seqnames %in% gsub("_","",isolate(input$chrom1)) &
    #                                                     anchorTwo(breakpoint_gint_full)@seqnames %in% gsub("_","",isolate(input$chrom2))]
    # 
    # if(
    #   grep(paste0("\\b",unique(as.character(submat_row_gr@seqnames)),"\\b"),gsub("_","",chromosomes))>grep(paste0("\\b",unique(as.character(submat_col_gr@seqnames)),"\\b"),gsub("_","",chromosomes))
    # ){
    #   SVs_data_in_submatrix_coords<-rebinGenomicInteractions(gint=breakpoint_gint_full_subset,
    #                                                          whole_genome_matrix = NULL,
    #                                                          rownames_gr = submat_col_gr,
    #                                                          colnames_gr = submat_row_gr,
    #                                                          rownames_mat = colnames(recast_matrix),
    #                                                          colnames_mat = rownames(recast_matrix),
    #                                                          method="nearest")
    # } else {SVs_data_in_submatrix_coords<-rebinGenomicInteractions(gint=breakpoint_gint_full_subset,
    #                                                                whole_genome_matrix = NULL,
    #                                                                rownames_gr = submat_row_gr,
    #                                                                colnames_gr = submat_col_gr,
    #                                                                rownames_mat = rownames(recast_matrix),
    #                                                                colnames_mat = colnames(recast_matrix),
    #                                                                method="nearest")
    # }
    #END CLIENT SIDE GINT PROCESSING
    # if(input$contour){
    #   p <- ggplot(data = ggplotmatrix, aes(x = Var2, y = Var1,fill=value,text=paste0("value:",value,"\nrow:",Var1,"\ncol:",Var2,"\n",value1)) ) + #geom_tile() + theme_void()
    #     geom_tile() + scale_x_discrete(breaks = block_index_labels_col) +
    #     scale_y_discrete(breaks = block_index_labels_row) + theme(axis.text.x = element_text(angle=60, hjust=1)) + 
    #     ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.5, limits = c(0, 1)) +  theme(legend.position="bottom",axis.title = element_blank()) + coord_flip() #+ scale_y_reverse(breaks=block_indices)
    #   }
    #rep(paste0(colnames(lumpy_summarized_counts[,3:ncol(lumpy_summarized_counts)]),collapse='/n'),nrow(lumpy_summarized_counts))
    #tidyr::unite(data = lumpy_summarized_counts[,3:ncol(lumpy_summarized_counts)],sep="\n")[,1]
    #
    if(exists("lumpy_summarized_counts") && isolate(input$lumpy_points_toggle)){
      
      lumpy_summarized_counts$textlabel<-unlist(strsplit(x = paste0("col:",lumpy_summarized_counts$row_bin_label,"\nrow:",lumpy_summarized_counts$col_bin_label,"\ntotal SVs:",lumpy_summarized_counts$total_samples,
                                                                    "\nhighest freq SV type:",lumpy_summarized_counts$highest_count_sample,lumpy_summarized_counts$highest_count_sample_count/lumpy_summarized_counts$total_samples*100,"%\n types, ranked:",lumpy_summarized_counts$concatenated_sample_names,collapse="@"),"@"))
      # p<-p + geom_point(data=lumpy_summarized_counts,mapping=aes(x=as.integer(as.character(lumpy_summarized_counts$col_bin_index)),y=as.integer(as.character(lumpy_summarized_counts$row_bin_index)),
      #                                                            color=lumpy_summarized_counts$highest_count_sample,size=lumpy_summarized_counts$total_samples,
      #                                                            text=lumpy_summarized_counts$textlabel
      # 
      #                                                            ))
      if(is.null(lumpy_summarized_counts$start1))
      {lumpy_summarized_counts<-dplyr::bind_cols(lumpy_summarized_counts,reshape2::colsplit(lumpy_summarized_counts$row_bin_label,"_",c("chr1","start1","end1")))
      lumpy_summarized_counts<-dplyr::bind_cols(lumpy_summarized_counts,reshape2::colsplit(lumpy_summarized_counts$col_bin_label,"_",c("chr2","start2","end2")))
      }
      p<-p + geom_point(data=lumpy_summarized_counts,mapping=aes(x=as.numeric(as.character(lumpy_summarized_counts$start1)),y=as.numeric(as.character(lumpy_summarized_counts$start2)),
                                                                 color=as.character(lumpy_summarized_counts$highest_count_sample),size=as.numeric(as.character(lumpy_summarized_counts$total_samples)),
                                                                 text=lumpy_summarized_counts$textlabel))
    }
    #
    if(exists("SVs_data_in_submatrix_coords") && isolate(input$plot_points_toggle))
    { SVs_data_in_submatrix_coords$col_bin_index<-as.numeric(as.character(SVs_data_in_submatrix_coords$col_bin_index))
    SVs_data_in_submatrix_coords$row_bin_index<-as.numeric(as.character(SVs_data_in_submatrix_coords$row_bin_index))
    if(is.null(SVs_data_in_submatrix_coords$start1))
    {SVs_data_in_submatrix_coords<-dplyr::bind_cols(SVs_data_in_submatrix_coords,reshape2::colsplit(SVs_data_in_submatrix_coords$row_bin_label,"_",c("chr1","start1","end1")))
    SVs_data_in_submatrix_coords<-dplyr::bind_cols(SVs_data_in_submatrix_coords,reshape2::colsplit(SVs_data_in_submatrix_coords$col_bin_label,"_",c("chr2","start2","end2")))
    }
    
    SVs_data_in_submatrix_coords$textlabel<-unlist(strsplit(x = paste0("col:",SVs_data_in_submatrix_coords$row_bin_label,"\nrow:",SVs_data_in_submatrix_coords$col_bin_label,"\ntotal SVs:",SVs_data_in_submatrix_coords$total_samples,
                                                                       "\nhighest freq SV type:",SVs_data_in_submatrix_coords$highest_count_sample,SVs_data_in_submatrix_coords$highest_count_sample_count/SVs_data_in_submatrix_coords$total_samples*100,"%\n types, ranked:",
                                                                       SVs_data_in_submatrix_coords$concatenated_sample_names,collapse="@"),"@"))
    
    #print(p_with_points)
    #},error = function(err) {
    #                                                                                    print(paste("Caught & handled error:  ",err))
    
    tryCatch( highest_over_tot<-as.numeric(SVs_data_in_submatrix_coords$highest_count_sample_count/SVs_data_in_submatrix_coords$total_samples),error = function(e) NULL) 
    tryCatch(colorvals<-as.character(cut(highest_over_tot,breaks=unique(quantile(highest_over_tot,probs=c(0.25,0.5,0.75))))),error = function(e) NULL)  
    
    if(exists("colorvals"))
    {    p_with_points<-p + geom_point(data=SVs_data_in_submatrix_coords,mapping = aes(x=as.numeric(as.character(SVs_data_in_submatrix_coords$start1)),y=as.numeric(as.character(SVs_data_in_submatrix_coords$start2)),
                                                                                       text=SVs_data_in_submatrix_coords$textlabel,
                                                                                       size=as.numeric(as.character(SVs_data_in_submatrix_coords$total_samples)),
                                                                                       #shape=as.character(SVs_data_in_submatrix_coords$highest_count_sample),
                                                                                       color= colorvals) ) + labs(color="",size="")
    } else {
      p_with_points<-p + geom_point(data=SVs_data_in_submatrix_coords,mapping = aes(x=as.numeric(as.character(SVs_data_in_submatrix_coords$start1)),y=as.numeric(as.character(SVs_data_in_submatrix_coords$start2)),
                                                                                    text=SVs_data_in_submatrix_coords$textlabel,
                                                                                    color="CGI SV",
                                                                                    size=as.numeric(as.character(SVs_data_in_submatrix_coords$total_samples))) ) + labs(size="")}
    
    
    
    #+ scale_color_gradient(low="green",high="darkgreen") 
    #color=as.numeric(SVs_data_in_submatrix_coords$highest_count_sample_count/SVs_data_in_submatrix_coords$total_samples)
    # + scale_colour_gradientn(colours = c("blue","white","red"),values=c(0,0.5,1)) 
    # p_with_points<-p + geom_point(data=SVs_data_in_submatrix_coords,mapping = aes(x=as.integer(as.character(SVs_data_in_submatrix_coords$col_bin_index)),y=as.integer(as.character(SVs_data_in_submatrix_coords$row_bin_index)),
    #                                                                               text=tidyr::unite(data = SVs_data_in_submatrix_coords[,3:ncol(SVs_data_in_submatrix_coords)],sep="\n")[,1],
    #                                                                               size=as.integer(as.character(SVs_data_in_submatrix_coords$total_samples)),
    #                                                                               #shape=as.character(SVs_data_in_submatrix_coords$highest_count_sample),
    #                                                                               color=as.character(arules::discretize(as.numeric(SVs_data_in_submatrix_coords$highest_count_sample_count/SVs_data_in_submatrix_coords$total_samples),method="interval"))
    #                                                                               #color=as.numeric(SVs_data_in_submatrix_coords$highest_count_sample_count/SVs_data_in_submatrix_coords$total_samples)
    # )) #+ scale_colour_gradientn(colours = c("blue","white","red"),values=c(0,0.5,1)) 
    
    #                                      scale_colour_gradient2()
    #set the range to be specific if there are coordinates (the cell +/- 4), else choose the max range for the particular axis.
   if(debug){browser()}
    
    
    #check for the correct format.
    plotly_output<-plotly::ggplotly(p_with_points,tooltip="text") %>% layout(margin=list(r=0, l=200, t=0, b=200),width=isolate(input$heatmapHeight),height=round(isolate(input$heatmapHeight)/1.25))
    } else {if(exists("p"))
    {
      plotly_output<-plotly::ggplotly(p,tooltip="text") %>% layout(margin=list(r=0, l=200, t=0, b=200),width=isolate(input$heatmapHeight),height=round(isolate(input$heatmapHeight)/1.25))
    }
    }
    #
    
    #plotly_output<-plotly::ggplotly(p) %>%       layout(margin=list(r=0, l=200, t=0, b=200),width=1280,height=1024)
    #%>% saveWidget(title = gsub("_","",paste0(chromosomes[isolate(input$chrom1)],"-",chromosomes[isolate(input$chrom2)])),file = paste0(chromosomes[isolate(input$chrom1)],chromosomes[isolate(input$chrom2)],"transparent_tooltipv27_coord_no_flip_downsample_upward_orientation_plotly_nrsample.html"),selfcontained = T)
    #
    if( (!is.null(isolate(input$loc_input_row)) | !is.null(isolate(input$loc_input_col)) ) & (!isolate(input$loc_input_row)=="" | !isolate(input$loc_input_col)==""))
    {
     if(debug){browser()}
      #acknowledgement: thanks to stackoverflow comments that made package a reality.
      #find the location of the bin in terms of map coordinates for x
      #store this as the xcentercoord
      #do the same for y
      #store as ycentercoord
      rowsplit<-reshape2::colsplit(isolate(input$loc_input_row),c("\\:|\\-"),c("chr","start","end"))
      columnsplit<-reshape2::colsplit(isolate(input$loc_input_col),c("\\:|\\-"),c("chr","start","end"))
      xmin<-columnsplit$start
      xmin<-ggplotmatrix$start2[which.min(abs(ggplotmatrix$start2-xmin))]-1e6
      xmax<-columnsplit$end
      xmax<-ggplotmatrix$start2[which.min(abs(ggplotmatrix$start2-xmax))]+1e6
      ymin<-rowsplit$start
      ymin<-ggplotmatrix$start2[which.min(abs(ggplotmatrix$start1-ymin))]-1e6
      ymax<-rowsplit$end
      ymax<-ggplotmatrix$start2[which.min(abs(ggplotmatrix$start1-ymax))]+1e6
      xglobalmin<-min(ggplotmatrix$start2)
      yglobalmin<-min(ggplotmatrix$start1)
      xglobalmax<-max(ggplotmatrix$start2)
      yglobalmax<-max(ggplotmatrix$start1)
      #edge case-- if the xcentercoord is greater than max or less than zero, reset to zero.
      #edge case-- do the same for y
      
      if(xmin<xglobalmin){xmin<-xglobalmin}
      if(ymin<yglobalmin){ymin<-yglobalmin}
      #edge case-- if xmax is greater than the maximum y, then reset to max.
      #edge case-- do the same for y
      
      if(xmax>xglobalmax){xmax<-xglobalmax}
      if(ymax>yglobalmax){ymax<-yglobalmax}
      #ggplotly(p, dynamicTicks = T) %>% layout(xaxis=list(autorange=F, range=c(xcentercoord-4,xcentercoord+4)), yaxis=list(autorange=F, range=c(20,30)))
      if(!exists("xmin")){xmin<-xglobalmin}
      if(!exists("xmax")){xmax<-xglobalmax}
      if(!exists("ymin")){ymin<-yglobalmin}
      if(!exists("ymax")){ymax<-yglobalmax} #need to round the max and min for all.
      #xmin<-floor(xmin/1e6)*1e6
      
      if(exists("p_with_points")){
        plotly_output<-plotly::ggplotly(p_with_points,tooltip="text") %>% layout(margin=list(r=0, l=200, t=0, b=200),width=isolate(input$heatmapHeight),height=round(isolate(input$heatmapHeight)/1.25),
                                                                                 xaxis=list(range=c(xmin,xmax),autorange=F), yaxis=list(range=c(ymin,ymax),autorange=F))
      } else {
        if(exists("p"))
        {

          plotly_output<-plotly::ggplotly(p,tooltip="text") %>% layout(margin=list(r=0, l=200, t=0, b=200),width=isolate(input$heatmapHeight),height=round(isolate(input$heatmapHeight)/1.25),xaxis=list(range=c(xmin,xmax),autorange=F), yaxis=list(range=c(ymin,ymax),autorange=F))
        }
      }
      return(plotly_output)
    } else {}
    
   if(debug){browser()}
    print(plotly_output)
  })
  outputOptions(output,"plotlyChromosomalHeatmap",suspendWhenHidden=F)
  output$whole_genome_image<-renderUI({
    
    input$whole_genome_max_cap
    input$goButton
    #browser()
    if(isolate(input$data_source)=="linreg_osteosarcoma_CNVkit")
    {
      data_prefix<-"osteo"
    }
    if(isolate(input$data_source)=="TCGA_NBL_low_pass")
    {
      data_prefix<-"nbl"
    }
    if(is.null(data_prefix)){return(NULL)}
    # list(src = paste0("http://alps.nci.nih.gov/james/plotly_dashboard/whole_genome_pngs/",data_prefix,"_whole_genome_full_no_downsample_no_labels_rescaled_max_cap_",isolate(input$whole_genome_max_cap),".png"),
    #      contentType = 'image/png',
    #      width = isolate(input$heatmapHeight),
    #      height = round(isolate(input$heatmapHeight)/1.25),
    #      alt = "This is alternate text")
    tags$img(src = paste0("http://alps.nci.nih.gov/james/plotly_dashboard/whole_genome_pngs/",data_prefix,"_whole_genome_full_no_downsample_no_labels_rescaled_max_cap_",isolate(input$whole_genome_max_cap),".png"),
             #           contentType = 'image/png',
             width = isolate(input$heatmapHeight),
             height = round(isolate(input$heatmapHeight)/1.25),
             alt = "whole genome png")
  }) 
  
  # output$freq_table<-renderDataTable({
  # 
  #  return(data.table()) 
  # })
  getGGplotMatrix<-function(){if(exists("ggplotmatrix")){return(ggplotmatrix)}else{return(NULL)}}
  getGGplotMatrix_full<-function(){if(exists("ggplotmatrix_full")){return(ggplotmatrix_full)}else{return(NULL)}}
  
  #TCGA_low_pass_sample_info
  get_tcga_lp_sample_info<-function(){if(exists("TCGA_low_pass_sample_info")){return(TCGA_low_pass_sample_info)}else{return(NULL)}}
  get_recast_matrix<-function(){if(exists("recast_matrix")){return(recast_matrix)}else{return(NULL)}}
  get_downsample_factor<-function(){if(exists("downsample_factor")){return(downsample_factor)}else{return(NULL)}}
  get_recast_matrix_full<-function(){if(exists("recast_matrix_full")){return(recast_matrix_full)}else{return(NULL)}}
  get_rownames_gr<-function(){if(exists("rownames_gr")){return(rownames_gr)}else{return(NULL)}}
  get_colnames_gr<-function(){if(exists("colnames_gr")){return(colnames_gr)}else{return(NULL)}}
  get_rownames_gr_full<-function(){if(exists("rownames_gr_full")){return(rownames_gr_full)}else{return(NULL)}}
  get_colnames_gr_full<-function(){if(exists("colnames_gr_full")){return(colnames_gr_full)}else{return(NULL)}}
  
  # get_recast_matrix<-function(){return(recast_matrix)}
  output$expression_data<-DT::renderDataTable({
    #browser()
    if(is.null(event_data("plotly_click"))){return(data.table())}
    if(isolate(input$data_source)=="linreg_osteosarcoma_CNVkit")
    {     
      
      recast_matrix<-get_recast_matrix()
      row_label<-rownames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1] #correct column label.
      column_label<-colnames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1] #correct column label.
      #row_point_gr<-underscored_pos_to_GRanges(row_label)
      #column_point_gr<-underscored_pos_to_GRanges(column_label)
      #row_index<-as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1
      #col_index<-as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1 #row and col indices of the subset matrix.
      row_index_full<-grep(row_label,rownames(get_recast_matrix_full()))
      col_index_full<-grep(column_label,colnames(get_recast_matrix_full()))
      #
      #rowclick<-length(common_coords)-myReactives$currentClick$lat
      #colclick<-myReactives$currentClick$lng
     if(debug){browser()}
      rowexpression<-as.data.table(subsetByOverlaps(expression_data_gr,get_rownames_gr_full()[seq(from=row_index_full,to=row_index_full+3)]))
      colexpression<-as.data.table(subsetByOverlaps(expression_data_gr,get_colnames_gr_full()[seq(from=col_index_full,to=col_index_full+3)]))} else {
        if(isolate(input$data_source)=="TCGA_NBL_low_pass" | isolate(input$data_source) %in% c("TCGA_NBL_stage3_subset","TCGA_NBL_stage4_subset","TCGA_NBL_stage4s_subset","TCGA_NBL_myc_amp_subset","TCGA_NBL_not_myc_amp_subset"))
        {
          
if(debug){browser()}
          rownames_gr_full<-get_rownames_gr_full()
          colnames_gr_full<-get_colnames_gr_full()
#         if(!exists("expression_data_gr_nbl")){
            tryCatch(expression_data_gr_nbl<-readRDS(paste0(get("basefn",.GlobalEnv),"tcga_nbl_expression.rds")),error = function(e) NULL)  
 #         }
          if(length(expression_data_gr_nbl)==0){
          tryCatch(expression_data_gr_nbl<-readRDS(paste0(get("basefn",.GlobalEnv),"tcga_nbl_expression.rds")),error = function(e) NULL)
            }
          #mcols(expression_data_gr_nbl)$SYMBOL<-expression_data_gr_nbl$....external_gene_name
         if(debug){browser()}
          rowexpression<-as.data.table(subsetByOverlaps(expression_data_gr_nbl,rownames_gr_full[rownames_gr_full@ranges@start==event_data("plotly_click")[["y"]]]))
          colexpression<-as.data.table(subsetByOverlaps(expression_data_gr_nbl,colnames_gr_full[colnames_gr_full@ranges@start==event_data("plotly_click")[["x"]]]))
        }
      }
    
    rowexpression$rowcol<-"row"
    colexpression$rowcol<-"col"
    comb_expression_df<-rbind(rowexpression,colexpression)
    #comb_expression_df_t<-as.data.table(t(comb_expression_df))
    #return(comb_expression_df_t)
    # cat(file=stderr(),paste0("expression_data"))
    # cat(file=stderr(),ls())
    #make the rownames match for nbl
    outputexpression_df<-as.data.table(unique(comb_expression_df[,c("SYMBOL","seqnames","start","end","gene_type","rowMean","rowMeanPctl","rowVar","rowVarPctl")]))
    outputexpression_df_sorted<-outputexpression_df[order(-outputexpression_df$rowVarPctl),]
    return(as.data.table(outputexpression_df_sorted))
  })
  output$census_data<-DT::renderDataTable({
    #
    if(is.null(event_data("plotly_click"))){return(data.table())}
    recast_matrix<-get_recast_matrix()
    if(length(intersect(ls(),"census_data_gr"))!=1) {    tryCatch(census_data_gr<-readRDS(paste0(basefn,"censushg19.rds")),error = function(e) NULL)}
    row_label<-rownames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1] #correct column label.
    column_label<-colnames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1] #correct column label.
    #row_point_gr<-underscored_pos_to_GRanges(row_label)
    #column_point_gr<-underscored_pos_to_GRanges(column_label)
    #row_index<-as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1
    #col_index<-as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1 #row and col indices of the subset matrix.
    row_index_full<-grep(row_label,rownames(get_recast_matrix_full()))
    col_index_full<-grep(column_label,colnames(get_recast_matrix_full()))
    #
    #rowclick<-length(common_coords)-myReactives$currentClick$lat
    #colclick<-myReactives$currentClick$lng
    rowcensus<-as.data.table(subsetByOverlaps(census_data_gr,get_rownames_gr_full()[seq(from=row_index_full,to=row_index_full+3)]))
    colcensus<-as.data.table(subsetByOverlaps(census_data_gr,get_colnames_gr_full()[seq(from=col_index_full,to=col_index_full+3)]))
    rowcensus$rowcol<-"row"
    colcensus$rowcol<-"col"
    comb_census_df<-rbind(rowcensus,colcensus)
    comb_census_df_t<-as.data.table(t(comb_census_df))
    # cat(file=stderr(),paste0("census_data"))
    # cat(file=stderr(),ls())
    #return(comb_census_df_t)
    #browser()
    return(unique(as.data.table(comb_census_df))) #[,c("SYMBOL","seqnames","start","end","gene_type","rowMean","rowMeanPctl","rowVar","rowVarPctl")]
  })
  
  # output$census_data<-renderDataTable({
  #   row_label<-rownames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1] #correct column label.
  #   column_label<-colnames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1] #correct column label.
  #   if(is.null(myReactives$currentClick)){return(data.frame())}
  #   #
  #   rowclick<-round(length(common_coords)-myReactives$currentClick$lat)
  #   colclick<-round(myReactives$currentClick$lng)
  #   rowcensus<-as.data.table(subsetByOverlaps(census_data_gr,rownames_gr[rowclick]))
  #   colcensus<-as.data.table(subsetByOverlaps(census_data_gr,colnames_gr[colclick]))
  #   rowcensus$rowcol<-"row"
  #   colcensus$rowcol<-"col"
  #   comb_expression_df<-rbind(rowcensus,colcensus)
  #   comb_expression_df_t<-t(comb_expression_df)
  #   return(comb_expression_df_t)
  #   
  # })
  output$gene_data <-
    renderPrint({
      if(is.null(event_data("plotly_click"))){return(data.table())}
      
      row_label<-rownames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1] #correct column label.
      column_label<-colnames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1] #correct column label.
      #if(myReactives)
      #
      #all_input<-isolate(input)
      # cat(file=stderr(),paste0("gene_data"))
      # cat(file=stderr(),ls())
      rowclick<-length(common_coords)-myReactives$currentClick$lat
      colclick<-myReactives$currentClick$lng
      row_genes<-genev[rowclick]
      col_genes<-genev[colclick]
      #
      output<-paste0("row genes:",as.character(genev[rowclick]),
                     "column genes:",as.character(genev[colclick]))
      return(output)
      
    })
  output$network <- renderVisNetwork({
    if (input$goButton == 0) {return()}
    
    input$goButton
    #browser()
    #     ggplotmatrix_filtered<-ggplotmatrix[ggplotmatrix$value > summary(heatmaply::percentize(ggplotmatrix$value))["3rd Qu."] | ggplotmatrix$value < summary(heatmaply::percentize(ggplotmatrix$value))["1st Qu."], ]
    #     ggplotmatrix_filtered<-ggplotmatrix[heatmaply::percentize(ggplotmatrix$value) > 0.9999 | heatmaply::percentize(ggplotmatrix$value) < 0.0001, ]
    ggplotmatrix_filtered<-ggplotmatrix_full[order(ggplotmatrix_full$value),]
    ggplotmatrix_filtered<-ggplotmatrix_filtered[c(1:(isolate(input$n_nodes)/2),(nrow(ggplotmatrix_filtered)-(isolate(input$n_nodes)/2)):nrow(ggplotmatrix_filtered)),]
    ggplotmatrix_filtered<-ggplotmatrix_filtered[as.character(ggplotmatrix_filtered$Var1)!=as.character(ggplotmatrix_filtered$Var2),]
    vertex.attrs<-list(name = unique(c(as.character(ggplotmatrix_filtered$Var1), as.character(ggplotmatrix_filtered$Var2))))
    
    edges<-rbind(as.character(ggplotmatrix_filtered$Var1),as.character(ggplotmatrix_filtered$Var2))
    weights<-ggplotmatrix_filtered$value
    G <- graph.empty(n = 0, directed = T)
    G <- add.vertices(G, length(vertex.attrs$name), attr = vertex.attrs)
    G <- add.edges(G, edges,weight=weights)
    G_connected<-delete.isolates(G)
    #    weights_discretized<-arules::discretize(E(G_connected)$weight)
    #     G_connected_D3<-networkD3::igraph_to_networkD3(G_connected,group = as.character(arules::discretize(strength(G_connected))))
    # forceNetwork(Links = G_connected_D3$links, Nodes = G_connected_D3$nodes, 
    #              Source = 'source', Target = 'target', 
    #              NodeID = 'name',Group='group',fontSize = 14,zoom=T)
    G_connected_vis<-visNetwork::toVisNetworkData(G_connected)
    G_connected_vis$edges$value<-G_connected_vis$edges$weight
    col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
    G_connected_vis$nodes$color<-sapply(col_fun(heatmaply::percentize(strength(G_connected)))  ,function(x) substr(x,start = 1,stop =  7))
    visNetwork::visNetwork(nodes = G_connected_vis$nodes,edges = G_connected_vis$edges,width = isolate(input$heatmapHeight),height = round(isolate(input$heatmapHeight)/1.25))  %>%
      visInteraction(hover = TRUE) %>%
      visEvents(hoverNode = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes);
                ;}")
})
  
  output$shiny_return <- DT::renderDataTable({
    input$current_node_id
    if(is.null(isolate(input$current_node_id))){return(data.table())}
    
    #browser()
    #DT::datatable(iris, options = list(lengthMenu = c(5, 30, 50), pageLength = 5)
    #paste0(input$current_node_id)
    return(as.data.table(ggplotmatrix[ggplotmatrix$Var1 %in% isolate(input$current_node_id) | ggplotmatrix$Var2 %in% isolate(input$current_node_id),]))#c("Var1,Var2","value","value1")
  },options = list(pageLength=5))#
  #pageLength = 5)
  output$sample_info<-renderPlotly({
    input$sample_hist_alpha
    if(is.null(event_data("plotly_click"))){return(data.table())}
    #browser()
    #ed <- event_data("plotly_click")
    if (is.null(event_data("plotly_click"))) {return("Click events appear here (double-click to clear)")}
    if(isolate(input$data_source)=="linreg_osteosarcoma_CNVkit" | isolate(input$data_source)=="TCGA_NBL_low_pass" | 
       isolate(input$data_source) %in% c("TCGA_NBL_stage3_subset","TCGA_NBL_stage4_subset","TCGA_NBL_stage4s_subset","TCGA_NBL_myc_amp_subset","TCGA_NBL_not_myc_amp_subset")
    )
    {
      recast_matrix<-get_recast_matrix()
      if(!is.null("recast_matrix")) {
        row_label<-rownames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1] #correct column label.
        column_label<-colnames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1] #correct column label.
        if(isolate(input$data_source)=="TCGA_NBL_low_pass" | 
           isolate(input$data_source) %in% c("TCGA_NBL_stage3_subset","TCGA_NBL_stage4_subset","TCGA_NBL_stage4s_subset","TCGA_NBL_myc_amp_subset","TCGA_NBL_not_myc_amp_subset"))
        {
          row_label<-paste0(isolate(input$chrom2),event_data("plotly_click")[["y"]],"_",event_data("plotly_click")[["y"]]+1e6-1)
          column_label<-paste0(isolate(input$chrom1),event_data("plotly_click")[["x"]],"_",event_data("plotly_click")[["x"]]+1e6-1)
        }
        if(length(bin_data$probe)==0)
        {
          bin_data$probe<-rownames(bin_data)
        }
        d<-as.data.table(bin_data[bin_data$probe %in% c(row_label,column_label),])
        if(nrow(d)==0){return("")}
        #p <- plotly::plot_ly(x = bin_data[1,], type = "histogram")
        #        cat(file=stderr(),paste0("sample_info"))
        #        cat(file=stderr(),ls())
        sample_info_p <- plot_ly(alpha = isolate(input$sample_hist_alpha)) %>%
          add_histogram(x = as.numeric(d[1,]),name=d[1,"probe"]) %>%
          add_histogram(x = as.numeric(d[2,]),name=d[2,"probe"]) %>%
          layout(barmode = "overlay")
        
        print(sample_info_p)
        return(sample_info_p)
      }
    } #end code for in-house data.
    if(isolate(input$data_source) %in% c("TCGA_AML_low_pass","TCGA_BRCA_low_pass","TCGA_OS_low_pass","TCGA_PRAD_low_pass"))
    {       
      TCGA_low_pass_sample_info<-get_tcga_lp_sample_info()
      
    }
  })
  output$sample_info_scatter<-renderPlotly({
    if(is.null(event_data("plotly_click"))){return(data.table())}
    #browser()
    if (is.null(event_data("plotly_click"))) {return("Click events appear here (double-click to clear)")}
    recast_matrix<-get_recast_matrix()
    if(!is.null("recast_matrix")) {
      row_label<-rownames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1] #correct column label.
      column_label<-colnames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1] #correct column label.
      if(isolate(input$data_source)=="TCGA_NBL_low_pass" | 
         isolate(input$data_source) %in% c("TCGA_NBL_stage3_subset","TCGA_NBL_stage4_subset","TCGA_NBL_stage4s_subset","TCGA_NBL_myc_amp_subset","TCGA_NBL_not_myc_amp_subset"))
      {
        row_label<-paste0(isolate(input$chrom2),event_data("plotly_click")[["y"]],"_",event_data("plotly_click")[["y"]]+1e6-1)
        column_label<-paste0(isolate(input$chrom1),event_data("plotly_click")[["x"]],"_",event_data("plotly_click")[["x"]]+1e6-1)
      }
      if(length(bin_data$probe)==0)
      {
        bin_data$probe<-rownames(bin_data)
      }
      
      d<-as.data.table(bin_data[bin_data$probe %in% c(row_label,column_label),])
      #testing
      #browser()
      bin_data_colsplit<-reshape2::colsplit(bin_data$probe,"_",c("chr","start","end"))
      bin_data_colsplit[bin_data_colsplit$chr=="chr19",]
      #end testing
      if(nrow(d)==0){return("")}
      #p <- plotly::plot_ly(x = bin_data[1,], type = "histogram")
      # cat(file=stderr(),paste0("census_data"))
      # cat(file=stderr(),ls())
      sample_info_p_scatter <- plot_ly(alpha = 0.6) %>%
        add_trace(x = as.numeric(d[1,]),name=d[1,"probe"],y=seq(1:ncol(d))) %>%
        add_trace(x = as.numeric(d[2,]),name=d[2,"probe"],y=seq(1:ncol(d)))# %>%
      # layout(barmode = "overlay")
      print(sample_info_p_scatter)
    }
    
  })
  output$minimap<-renderPlotly({
    if(is.null(event_data("plotly_click"))){return(data.table())}
    
    if (is.null(event_data("plotly_click"))) {return("Click events appear here (double-click to clear)")}
    recast_matrix<-get_recast_matrix()
    ggplotmatrix_full<-getGGplotMatrix_full()
    recast_matrix_full<-get_recast_matrix_full()
    if(!is.null("recast_matrix") & !is.null("recast_matrix_full")) {
      row_label<-rownames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1] #correct column label.
      column_label<-colnames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1] #correct column label.
      if(isolate(input$data_source)=="TCGA_NBL_low_pass" | 
         isolate(input$data_source) %in% c("TCGA_NBL_stage3_subset","TCGA_NBL_stage4_subset","TCGA_NBL_stage4s_subset","TCGA_NBL_myc_amp_subset","TCGA_NBL_not_myc_amp_subset"))
      {
        row_label<-paste0(isolate(input$chrom2),event_data("plotly_click")[["y"]],"_",event_data("plotly_click")[["y"]]+1e6-1)
        column_label<-paste0(isolate(input$chrom1),event_data("plotly_click")[["x"]],"_",event_data("plotly_click")[["x"]]+1e6-1)
      }
      if(length(bin_data$probe)==0)
      {
        bin_data$probe<-rownames(bin_data)
      }
      d<-as.data.table(bin_data[bin_data$probe %in% c(row_label,column_label),])
      if(nrow(d)==0){return("")}
      row_labels_minimap<-rownames(recast_matrix_full)[grep(row_label,rownames(recast_matrix_full)):(grep(row_label,rownames(recast_matrix_full))+3)] #we subset by every fourth number along the rows and columns, hence we need n, n+1, n+2, n+3 (or n1:n2-1, the first number and all the numbers leading up to the next).
      col_labels_minimap<-colnames(recast_matrix_full)[grep(column_label,colnames(recast_matrix_full)):(grep(column_label,colnames(recast_matrix_full))+3)]
      ggplotmatrix_minimap<-ggplotmatrix_full[as.character(ggplotmatrix_full$Var1) %in% row_labels_minimap & as.character(ggplotmatrix_full$Var2) %in% col_labels_minimap, ]
      p <- ggplot(data = ggplotmatrix_minimap ) + #geom_tile() + theme_void()
        geom_raster(aes(x = Var2, y = Var1,fill=value,text=paste0("value:",value,"\nrow:",Var1,"\ncol:",Var2,"\n",value1))) + scale_x_discrete() +
        scale_y_discrete() + theme(axis.text.x = element_text(angle=60, hjust=1)) + 
        ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.5, limits = c(0, 1)) +  theme(legend.position="bottom",axis.title = element_blank()) #+ coord_flip() #+ scale_y_reverse(breaks=block_indices)
      # cat(file=stderr(),paste0("minimap"))
      # cat(file=stderr(),ls())
      
      plotly_output<-plotly::ggplotly(p,tooltip="text") %>% layout(margin=list(r=0, l=200, t=0, b=200),width=isolate(input$heatmapHeight),height=isolate(input$heatmapHeight)/1.25)
      #print(plotly_output)
      #essentially, grab the row and column bins (above) for the sampled matrix, then grab the same coordinates for the full matrix, plus four to x, plus four to y.
      #p <- plotly::plot_ly(x = bin_data[1,], type = "histogram")
      # sample_info_p_scatter <- plot_ly(alpha = 0.6) %>%
      #   add_trace(x = as.numeric(d[1,]),name=d[1,"probe"],y=seq(1:ncol(d))) %>%
      #   add_trace(x = as.numeric(d[2,]),name=d[2,"probe"],y=seq(1:ncol(d)))# %>%
      # # layout(barmode = "overlay")
      # print(sample_info_p_scatter)
    }
  })
  output$sample_info_scatter2<-renderPlotly({
    
    if (is.null(event_data("plotly_click"))) {return(data.table())}
    #browser()
    if(isolate(input$data_source)=="linreg_osteosarcoma_CNVkit" | isolate(input$data_source)=="TCGA_NBL_low_pass" | 
       isolate(input$data_source) %in% c("TCGA_NBL_stage3_subset","TCGA_NBL_stage4_subset","TCGA_NBL_stage4s_subset","TCGA_NBL_myc_amp_subset","TCGA_NBL_not_myc_amp_subset"))
    {
      recast_matrix<-get_recast_matrix()
      if(!is.null("recast_matrix")) {
        #
        row_label<-rownames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1] #correct column label.
        column_label<-colnames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1] #correct column label.
        if(isolate(input$data_source)=="TCGA_NBL_low_pass" | 
           isolate(input$data_source) %in% c("TCGA_NBL_stage3_subset","TCGA_NBL_stage4_subset","TCGA_NBL_stage4s_subset","TCGA_NBL_myc_amp_subset","TCGA_NBL_not_myc_amp_subset"))
        {
          row_label<-paste0(isolate(input$chrom2),event_data("plotly_click")[["y"]],"_",event_data("plotly_click")[["y"]]+1e6-1)
          column_label<-paste0(isolate(input$chrom1),event_data("plotly_click")[["x"]],"_",event_data("plotly_click")[["x"]]+1e6-1)
        }
        if(length(bin_data$probe)==0)
        {
          bin_data$probe<-rownames(bin_data)
        }
        
        d<-as.data.table(bin_data[bin_data$probe %in% c(row_label,column_label),])
        if(nrow(d)==0){return("")}
        d_sample_names<-names(d)[2:length(names(d))]
        #p <- plotly::plot_ly(x = bin_data[1,], type = "histogram")
        #
        #names(d)
        #sample_info_p_scatter2 <- plot_ly(alpha = 0.6,x = as.numeric(d[1,]),y=as.numeric(d[2,]),name=d_sample_names)
        #
        
        d_t<-as.data.frame(t(d)[2:ncol(d),])
        colnames(d_t)<-d$probe
        d_t<-as.data.frame(sapply(as.data.frame(d_t),function(x) as.numeric(as.character(x))))
        rownames(d_t)<-d_sample_names
        if(ncol(d_t)==1){d_t[,2]<-d_t[,1]
        colnames(d_t)[2]<-paste0(d$probe,"_")}
        #,text=paste0("x: ",paste0(colnames(d_t)[1]),"  ", d_t[,1],"\n y:",paste0(colnames(d_t)[2]),"  ",d_t[,2],"\n ",rownames(d_t))
        #,color=rownames(d_t)
        #
        sample_info_p_scatter2<-ggplot(data = d_t,aes(x=d_t[,1],y=d_t[,2])) + geom_point(aes(color=rownames(d_t),text=paste0("x: ",paste0(colnames(d_t)[1]),"  ", d_t[,1],"\n y:",paste0(colnames(d_t)[2]),"  ",d_t[,2],"\n ",rownames(d_t)))) + theme(legend.position="none") +
          xlab(paste0(colnames(d_t)[1])) + ylab(paste0(colnames(d_t)[2])) + geom_smooth(method=lm)
        # %>% #name=d[1,"probe"],y=seq(1:ncol(d))
        #add_trace(x = as.numeric(d[2,]),name=d[2,"probe"],y=seq(1:ncol(d)))# %>%
        # layout(barmode = "overlay")
        
        # cat(file=stderr(),paste0("sample_info_scatter2"))
        #cat(file=stderr(),ls())
        #     cat(file=stderr(),sapply(ls(),function(x) paste0(unlist(paste0(head(get(x)))))))
        # cat(file=stderr(),paste0("sample_info_p_scatter2"))
        # cat(file=stderr(),str(sample_info_p_scatter2))
        # cat(file=stderr(),paste0("sample_info_p_scatter2_length"))
        # cat(file=stderr(),length(sample_info_p_scatter2))
        # cat(file=stderr(),unlist(sapply(ls(),function(x) paste0(paste0(head(get(x)))))))
        # cat(file=stderr(),paste0("sample_info_p_scatter2"))
        # cat(file=stderr(),str(sample_info_p_scatter2))
        
        
        #cat(file=stderr(),sapply(ls(),function(x) get(x)))
        print(plotly::ggplotly(sample_info_p_scatter2,tooltip=c("text")))
        return(plotly::ggplotly(sample_info_p_scatter2,tooltip=c("text")))
      }
    } #end in-house data processing
    
    if(isolate(input$data_source) %in% c("TCGA_AML_low_pass","TCGA_BRCA_low_pass","TCGA_OS_low_pass","TCGA_PRAD_low_pass"))
    {
      
      TCGA_low_pass_sample_info<-get_tcga_lp_sample_info()
      recast_matrix <- get_recast_matrix()
      if (!is.null("recast_matrix")) {
        row_label <- rownames(recast_matrix)[order(get_rownames_gr())][as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1])) + 1]
        column_label <- colnames(recast_matrix)[order(get_colnames_gr())][as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2])) + 1]
        if(isolate(input$data_source)=="TCGA_NBL_low_pass")
        {
          row_label<-paste0(isolate(input$chrom2),event_data("plotly_click")[["y"]],"_",event_data("plotly_click")[["y"]]+1e6-1)
          column_label<-paste0(isolate(input$chrom1),event_data("plotly_click")[["x"]],"_",event_data("plotly_click")[["x"]]+1e6-1)
        }
        d<-as.data.table(TCGA_low_pass_sample_info[TCGA_low_pass_sample_info$pos %in% c(row_label,column_label),])
        if("TCGA_CNV_data_gr.....relativeCvg" %in% colnames(TCGA_low_pass_sample_info)){
          d<-as.data.table(TCGA_low_pass_sample_info[TCGA_low_pass_sample_info$pos %in% c(row_label,column_label),c("TCGA_CNV_data_gr.....relativeCvg","TCGA_CNV_data_gr.....sample")])
          d_row<-as.data.table(TCGA_low_pass_sample_info[TCGA_low_pass_sample_info$pos %in% c(row_label),c("TCGA_CNV_data_gr.....relativeCvg","TCGA_CNV_data_gr.....sample")])
          d_col<-as.data.table(TCGA_low_pass_sample_info[TCGA_low_pass_sample_info$pos %in% c(column_label),c("TCGA_CNV_data_gr.....relativeCvg","TCGA_CNV_data_gr.....sample")])
        } else { if("TCGA_CNV_data_gr.....relativeCvg" %in% colnames(TCGA_low_pass_sample_info))
          d<-as.data.table(TCGA_low_pass_sample_info[TCGA_low_pass_sample_info$pos %in% c(row_label,column_label),c("TCGA_CNV_data_gr.....Segment_Mean","TCGA_CNV_data_gr.....sample")])
        d_row<-as.data.table(TCGA_low_pass_sample_info[TCGA_low_pass_sample_info$pos %in% c(row_label),c("TCGA_CNV_data_gr.....Segment_Mean","TCGA_CNV_data_gr.....sample")])
        d_col<-as.data.table(TCGA_low_pass_sample_info[TCGA_low_pass_sample_info$pos %in% c(column_label),c("TCGA_CNV_data_gr.....Segment_Mean","TCGA_CNV_data_gr.....sample")])
        }
        
        if(nrow(d)==0){return("")}
        sample_info_p_scatter2<-ggplot(data = d_row,aes(x=unlist(d_row[,1]),y=unlist(d_col[,1]))) + 
          geom_point(aes(color=unlist(d_row[,2]),shape=unlist(d_col[,2]),
                         text=paste0("row_value: ",paste0(d_row[,1]),"/n sample: ",paste0(d_row[,2]),
                                     " col_value: ", d_col[,1],"\n sample:",paste0(d_col[,2])))) + theme(legend.position="none") +
          xlab("column segmentation value") + ylab("row segmentation value") + geom_smooth(method=lm)
        # cat(file=stderr(),paste0("sample_info_scatter2"))
        # cat(file=stderr(),ls())
        
      }
      #d["TCGA_CNV_data_gr.....sample"
    }
  })
  output$freq_table <- DT::renderDataTable({
    #if(isolate(is.null(input$subset))){selected_rows<-1:nrow(mappability_df)} 
    #textv_subset<-textv[selected_rows]
    #d<-as.character(names(event_data("plotly_hover")))
    #
    # cat(file=stderr(),paste0(get_recast_matrix()
    #                          [
    #                            as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1,
    #                            as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1
    #                          ]))
    # cat(file=stderr(),rownames(get_recast_matrix())[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1])
    recast_matrix<-get_recast_matrix()
    #cat(file=stderr(),paste0(d))
    if(!is.null("recast_matrix")) {
      row_label<-rownames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1] #correct column label.
      column_label<-colnames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1] #correct column label.
      d<-as.data.table(freq_data[freq_data$pos %in% c(row_label,column_label)])
      # cat(file=stderr(),paste0("freq_table"))
      # cat(file=stderr(),ls())
      if (is.null(d)) {return(data.table())} else {
        return(d)}
    } else {return(data.table())}
    # cat(file=stderr(),paste0(event_data("plotly_click")))
    # cat(file=stderr(),paste0(names(event_data("plotly_click"))))
    # cat(file=stderr(),paste0(names(event_data("plotly_click")[["pointNumber"]])))
    # cat(file=stderr(),paste0(event_data("plotly_click")[["pointNumber"]]))
    # cat(file=stderr(),paste0(event_data("plotly_click")["pointNumber"]))
    # cat(file=stderr(),paste0(event_data("plotly_click")["curveNumber"]))
    # cat(file=stderr(),paste0(event_data("plotly_click")["x"]))
    #cat(file=stderr(),as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1 ) #row number
    #cat(file=stderr(),as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1 ) #col number
    # 
    # cat(file=stderr(),as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1 ) #col number
    # cat(file=stderr(),paste0(chromstarts_linreg))
    # cat(file=stderr(),paste0(head(common_coords_linreg)))
    # cat(file=stderr(),paste0(head(common_coords_linreg)))
    # cat(file=stderr(),paste(names(input)))
    # cat(file=stderr(),paste(input$chrom2))
    # cat(file=stderr(),paste(chromstarts_linreg[grep(input$chrom2,chromosomes)]+as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))))
    # cat(file=stderr(),paste(common_coords_linreg[chromstarts_linreg[grep(input$chrom2,chromosomes)]+as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))]))
    # cat(file=stderr(),paste(d))
    #need to convert to global coordinates
    #
    #cat(file=stderr(),exists("ggplotmatrix"))
    #cat(file=stderr(),exists("event_data(\"plotly_click\")"))
    #cat(file=stderr(),exists("event_data"))
    #cat(file=stderr(),paste0(event_data))
    #cat(file=stderr(),length(event_data))
    #cat(file=stderr(),paste0(event_data[[1]]))
    #cat(file=stderr(),paste0(signedRescale))
    # if(exists("ggplotmatrix") & !is.null(ggplotmatrix)){
    # recast_matrix<-reshape2::dcast(data=ggplotmatrix,formula=Var1 ~ Var2, var = ggplotmatrix$value) #this creates a matrix with 
    # if(ncol(recast_matrix)!=nrow(recast_matrix))
    # {
    #   rownames(recast_matrix)<-recast_matrix$Var1
    #   recast_matrix<-recast_matrix[,2:ncol(recast_matrix)]
    # }}
    # cat(file=stderr(),rownames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1] )
    # cat(file=stderr(),colnames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1] ) 
    # cat(file=stderr(),colnames(recast_matrix))
    # cat(file=stderr(),rownames(recast_matrix))
    # cat(file=stderr(),paste(head(ggplotmatrix)))
    # cat(file=stderr(),paste(input))
    # cat(file=stderr(),paste(names(input)))
    # cat(file=stderr(),paste0(chromosomes[as.integer(gsub("_","",gsub("chr","",isolate(input$chrom1))))]))
    
    # d<-freq_data[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1,as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1]
    #print(event_data("plotly_click"))
    #showLog()
    #class(event_data$plotly_click$pointNumber)
    #print(str(event_data("plotly_click")))
    #d<-as.data.table(event_data("plotly_click"))
    #d <-freq_data[as.integer(event_data("plotly_click")[["pointNumber"]]+1),]
    # if (is.null(d)) {return(data.table())} else {
    #   row_label<-rownames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][1]))+1]
    #   column_label<-colnames(recast_matrix)[as.integer(paste0(event_data("plotly_click")[["pointNumber"]][[1]][2]))+1]
    #   d<-as.data.table(freq_data[freq_data$pos %in% c(row_label,column_label)])
    # }
    # cat(file=stderr(),paste0(d))
    # return(d)
  })
}
