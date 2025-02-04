##' Start the Retistruct GUI
##' @return Object with \code{getData()} method to return
##' reconstructed retina data and environment \code{this} which
##' contains variables in object.
##' @importFrom shiny shinyApp
##' @export
retistruct <- function() {
  options(rgl.useNULL = TRUE) ## Prevents rgl from making its own window
  shinyApp(ui = ui, server = server)
}
