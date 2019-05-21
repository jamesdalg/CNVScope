#' Runs the CNVScope plotly shiny application.
#'
#' Runs the interactive suite of tools locally.
#' @name runCNVScopeLocal
#' @keywords CNV heatmap shiny plotly
#' @return none. Runs the application if the correct files are present.
#' @examples
#' \dontrun{
#' CNVScope::runCNVScopeLocal()
#' }
#' @export
runCNVScopeLocal<-function(){
  CNVScope::runCNVScopeShiny(useCNVScopePublicData = T)
}