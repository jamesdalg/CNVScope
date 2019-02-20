#' Runs the CNVScope plotly shiny application.
#'
#' Runs the interactive suite of tools locally or on a server if called in a script file (e.g. App.R).
#' @name runCNVScopeShiny
#' @keywords CNV heatmap shiny plotly
#' @import shinycssloaders shinythemes visNetwork ggplot2 reshape2 magrittr htmltools htmlwidgets jointseg logging foreach GenomicInteractions shinythemes
#' @importFrom BiocManager repositories
#' @rawNamespace import(circlize, except = degree)
#' @rawNamespace import(shiny, except = runExample)
#' @rawNamespace import(shinyjs, except = runExample)
#' @rawNamespace import(RCurl, except = reset)
#' @rawNamespace import(plotly, except = c(last_plot,select,filter))
#' @rawNamespace import(igraph, except = c(decompose, spectrum, groups))
#' @rawNamespace import(data.table, except = c(melt, dcast))
#' @rawNamespace import(GenomicFeatures ,except = show)
#' @param baseurl the url of the source files for the application (e.g. the contents of plotly_dashboard_ext). This will be pulled from remotely.
#' @param basefn the linux file path of the same source files.
#' @return none. Runs the application if the correct files are present.
#' @examples
#' \dontrun{
#' runCNVScopeShiny()
#' }
#' @export
globalVariables(c("common_coords_linreg","expression_data_gr","."), add=F)

runCNVScopeShiny<-function(baseurl=NULL,basefn=NULL) {
chrom.pairs<-NULL
  
options(scipen=999)
#if(getRversion() >= "2.15.1")  utils::globalVariables(c("."),add=F)
  
head.matrix<-function(mat,n=6L)
{
  return(mat[1:n,1:n])
}
delete.isolates <- function(graph, mode = 'all') {
  isolates <- which(igraph::degree(graph, mode = mode) == 0) 
  delete.vertices(graph, isolates)
}
options(repos = BiocManager::repositories())
options(shiny.error = browser)
options(shiny.fullstacktrace = TRUE)
getOption("repos")
options(shiny.sanitize.errors = F)
if(Sys.info()["nodename"]=="ncias-d2037-v.nci.nih.gov" | Sys.info()["nodename"]=="plotly.nci.nih.gov")
{
  baseurl<-"file:///srv/shiny-server/plotly_dashboard/"
  basefn<-"/srv/shiny-server/plotly_dashboard/"
} else {
  # baseurl<-"ftp://helix.nih.gov/pub/dalgleishjl/"
}
if(Sys.info()["nodename"]=="NCI-02105037-L")
{
  #baseurl<-"file:///W|/dalgleishjl/hicnv/"
  #baseurl<-"ftp://helix.nih.gov/pub/dalgleishjl/"
  baseurl<-"http://alps.nci.nih.gov/james/"
  #baseurl<-"file:///W:/dalgleishjl/hicnv/"
  basefn<-"W:/dalgleishjl/hicnv/"
}
baseurl<<-baseurl
basefn<<-basefn
tryCatch(bin_data<-readRDS((url(paste0(baseurl,"plotly_dashboard_ext/bin_data.rds")))),error = function(e) NULL) 
tryCatch(bin_data<-readRDS((paste0(basefn,"plotly_dashboard_ext/bin_data.rds"))),error = function(e) NULL) 
chromosomes<<-paste0("chr",c(seq(1:22),"X"),"_")
options(shiny.error = function() { 
  logging::logerror(sys.calls() %>% as.character %>% paste(collapse = ", ")) })
tryCatch(load(url(paste0(paste0(baseurl,"plotly_dashboard_ext/common_coords_linreg.RData")))),error = function(e) NULL)
tryCatch(load(paste0(paste0(basefn,"plotly_dashboard_ext/common_coords_linreg.RData"))),error = function(e) NULL)

chromstarts_linreg<-unlist(foreach(i=1:length(1:length(chromosomes))) %do% {grep(chromosomes[i],common_coords_linreg)[1]})
swap_row_col_genes=F
chrom.pairs<<-expand.grid(1:length(chromosomes),1:length(chromosomes))
chromosomes<-paste0("chr",c(seq(1:22),"X"),"_")
if(exists("basefn")) {#local objects:
  tryCatch(freq_data<-data.table::fread(paste0(basefn,"plotly_dashboard_ext/OS_freq_data.txt")),error = function(e) NULL)
  tryCatch(breakpoint_gint_full<-readRDS(paste0(basefn,"plotly_dashboard_ext/breakpoint_gint_full.rds")),error = function(e) NULL)
  tryCatch(expression_data_gr<-readRDS(paste0(basefn,"plotly_dashboard_ext/expression_data_gr.rds")),error = function(e) NULL)
  tryCatch(expression_data_gr_nbl<-readRDS(paste0(basefn,"plotly_dashboard_ext/tcga_nbl_expression.rds")),error = function(e) NULL)
  tryCatch(bin_data_gr<-readRDS(paste0(basefn,"plotly_dashboard_ext/bin_data_gr.rds")),error = function(e) NULL)
  #tryCatch(census_data_gr<-readRDS(paste0(basefn,"plotly_dashboard_ext/census_data_gr.rds")),error = function(e) NULL)
  tryCatch(census_data_gr<-readRDS(paste0(basefn,"plotly_dashboard_ext/censushg19.rds")),error = function(e) NULL)
  tryCatch(ensembl_gene_tx_data_gr<-readRDS(paste0(basefn,"plotly_dashboard_ext/ensembl_gene_tx_table_gr.rds")),error = function(e) NULL)
} else {
  if(exists("baseurl"))
  {tryCatch(freq_data<-data.table::fread(paste0(baseurl,"plotly_dashboard_ext/OS_freq_data.txt")),error = function(e) NULL)
    tryCatch(breakpoint_gint_full<-readRDS(url(paste0(baseurl,"plotly_dashboard_ext/breakpoint_gint_full.rds"))),error = function(e) NULL)
    tryCatch(expression_data_gr<-readRDS(url(paste0(baseurl,"plotly_dashboard_ext/expression_data_gr.rds"))),error = function(e) NULL)
    tryCatch(expression_data_gr_nbl<-readRDS(url(paste0(baseurl,"plotly_dashboard_ext/tcga_nbl_expression.rds"))),error = function(e) NULL)
    tryCatch(bin_data_gr<-readRDS(url(paste0(baseurl,"plotly_dashboard_ext/bin_data_gr.rds"))),error = function(e) NULL)
    #tryCatch(census_data_gr<-readRDS(url(paste0(baseurl,"plotly_dashboard_ext/census_data_gr.rds"))),error = function(e) NULL)
    tryCatch(census_data_gr<-readRDS(paste0(basefn,"plotly_dashboard_ext/censushg19.rds")),error = function(e) NULL)
    tryCatch(ensembl_gene_tx_data_gr<-readRDS(url(paste0(baseurl,"plotly_dashboard_ext/ensembl_gene_tx_table_gr.rds"))),error = function(e) NULL)
  }
}
CNVScopeui<-fluidPage(theme=shinytheme("flatly"), #shinythemes::themeSelector() 
                   tags$style(type="text/css",
                              ".shiny-output-error { visibility: hidden; }",
                              ".shiny-output-error:before { visibility: hidden; }"),
                   # Application title
                   titlePanel("Plotly Interchromosomal Heatmaps"),
                   
                   # Sidebar with a slider input for number of bins 
                   tabsetPanel(id = "tabs",tabPanel("Controls",fluidRow(column(width=2,offset = 1,
                                                                               #sidebarPanel(position="right",
                                                                               selectInput('data_source', 'data source', c("linreg_osteosarcoma_CNVkit","TCGA_NBL_low_pass","TCGA_NBL_stage3_subset","TCGA_NBL_stage4_subset","TCGA_NBL_stage4s_subset","TCGA_NBL_myc_amp_subset","TCGA_NBL_not_myc_amp_subset"), selected = "TCGA_NBL_low_pass"), #"TCGA_SARC_SNP6","TCGA_AML_low_pass","TCGA_BRCA_low_pass","TCGA_OS_low_pass" ,"TCGA_PRAD_low_pass"
                                                                               selectInput('chrom2', 'Chromosome (rows)', chromosomes, selected = "chr17_"),
                                                                               selectInput('chrom1', 'Chromosome (columns)', chromosomes, selected = "chr19_"),
                                                                               sliderInput('max_cap',"saturation limit",min=0.1,max=300,value = 75),
                                                                               #input$sample_hist_alpha
                                                                               sliderInput('heatmapHeight',"heatmap height",min=640,max=2048,value = 1024),
                                                                               #sliderInput('n_nodes',"number of nodes (top/bottom)",min=5,max=200,value = 50),
                                                                               conditionalPanel("input.data_source== 'linreg_osteosarcoma_CNVkit'",
                                                                                                checkboxInput('plot_points_toggle',"Plot Structural Variants",value = FALSE), 
                                                                                                checkboxInput('lumpy_points_toggle',"Plot Lumpy SVs",value = FALSE)),
                                                                               textInput('gene_input_row',"row_gene",NULL),
                                                                               textInput('loc_input_row',"row_location",NULL),
                                                                               textInput('gene_input_col',"col_gene",NULL),
                                                                               textInput('loc_input_col',"col_location",NULL),
                                                                               actionButton("geneSearch", "find genes"),
                                                                               actionButton("goButton", "create plots")
                                                                               
                                                                               
                   ))),tabPanel("Plots",fluidRow(column(5,offset=2,
                                                        h2("interactive chromosomal heatmap and minimap"),
                                                        withSpinner(plotlyOutput("plotlyChromosomalHeatmap",height = "100%"))#,
                                                        #               visNetwork::visNetworkOutput("network",height="1024")
                                                        
                                                        
                   ),column(1,offset=1,plotlyOutput("minimap",height=1024)#,verbatimTextOutput("shiny_return")
                   ))
                   ) ,#tabPanel("Top Network interactions", h2("Interactive Chromosomal Interaction network for top positive and negative relationships"),
                   #         fluidRow(dataTableOutput("shiny_return"),fluidRow(visNetwork::visNetworkOutput("network",height="1024"))#paste0(round(isolate(input$heatmapHeight)/1.25),"px")) 
                   #         )),
                   tabPanel("gain/loss frequency",
                            conditionalPanel("!is.null(event_data('plotly_click')) & is.null(output$freq_table)", fluidRow(h2("gain/loss frequency"), #
                                                                                                                           withSpinner(dataTableOutput("freq_table"))))),
                   tabPanel("COSMIC cancer gene census",h2("Cancer Gene Census Data"),
                            fluidRow( withSpinner(dataTableOutput("census_data")))), #end tabpanel
                   tabPanel("sample info",
                            fluidRow(column(2,offset=1,h3("sample histogram for row and column values at clicked point"),sliderInput('sample_hist_alpha',"histogram opacity",min=0.1,max=1,value = 0.6), withSpinner(plotlyOutput("sample_info"))),
                                     column(2,offset=1,h3("sample scatterplot for row and column segmentation values at clicked point"),withSpinner(plotlyOutput("sample_info_scatter"))),
                                     column(2,offset=1,h3("sample regression scatterplot for values at clicked point, colored by sample to show clustering"),withSpinner(plotlyOutput("sample_info_scatter2"))))      
                   ),
                   tabPanel("expression_data",
                            h2("expression data table for clicked point"),
                            fluidRow( withSpinner(dataTableOutput("expression_data")))
                   ),
                   tabPanel("Whole Genome View",fluidRow(column(11,offset=2,conditionalPanel("input.data_source== 'linreg_osteosarcoma_CNVkit' | input.data_source=='TCGA_NBL_low_pass'",h2("whole genome view"),
                                                                                             sliderInput(inputId = "whole_genome_max_cap",label = "Whole Genome p-value Saturation Cap",value = 75, min = 5,max=75,step = 5),
                                                                                             withSpinner(htmlOutput("whole_genome_image")))
                                                                
                   )#end column
                   )#end fluidRow
                   ) #end tabPanel for Whole Genome view
                   
                   ) #end tabset panel
                   
                   # Show a plot of the generated distribution
                   
                   
                   
                   
                   
)
shinyApp(ui = CNVScopeui, server = CNVScopeserver)
}
