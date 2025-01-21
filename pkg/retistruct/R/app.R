library(shiny)

# source("/R")
# source("../R/server.R")

##' Start the Retistruct GUI
##' @seealso gWidgets2
##' @return Object with \code{getData()} method to return
##' reconstructed retina data and environment \code{this} which
##' contains variables in object.
##' @importFrom grDevices dev.cur dev.set dev.print
##' @importFrom graphics identify
##' @importFrom utils packageVersion
##' @export
retistruct <- function() {
  options(rgl.useNULL = TRUE) ## Prevents rgl from making it's own window 
  shinyApp(ui = ui, server = server)
}
