##' Class to allow ggraphics objects to be magnified and scrolled 
##'
##' @description This creates a wrapper round
##'   \code{\link[gWidgets2]{ggraphics}} that adds + and - buttons to
##'   allow the user to resize the device. If the device is bigger
##'   than the plot window, then scrollbars appear.
##' @examples
##' \dontrun{
##' w <- gWidgets2::gwindow("Magnifier", height=400, width=400)
##' # Create the mangfier instance
##' m  <- Magnifier$new(container=w, ps=11)
##' ## Add extra buttons
##' g.print     <- gWidgets2::gbutton("Bitmap", container=m$buttons)
##' ## Plot in the window
##' plot(1:10, 1:10)
##' ## Set the device to be the magnifier
##' m$devSet()
##' plot(1:5, -5:-1)
##' ## Access the device with m$d
##' dev.set(m$d)
##' }
##' @author David Sterratt
Magnifier <- R6::R6Class(
  "Magnfier",
  public = list(
    ##' @field d Device handle
    d = NULL,
    ##' @field buttons gWidgets2::ggroup handler providing space for buttons
    buttons = NULL,
    ##' @description Magnifier constructor
    ##' @param container gWidgets2 parent container
    ##' @param width The width of the Magnifier in pixels
    ##' @param ... Arguments passed to \code{\link[gWidgets2]{ggraphics}}
    initialize = function(container, width=500, ...) {
      v1 <- gWidgets2::gvbox(cont=container, expand=TRUE)
      h  <- gWidgets2::ggroup(cont=v1)
      addSpace(h, width)
      self$buttons <- gWidgets2::ggroup(cont=v1)
      mag_plus <- gWidgets2::gbutton("+", cont=self$buttons)
      mag_minus <- gWidgets2::gbutton("-", cont=self$buttons)
      v2 <- gWidgets2::ggroup(horizontal=FALSE, cont=v1, use.scrollwindow=TRUE, expand=TRUE)
      gg <- gWidgets2::ggraphics(cont=v2, expand=TRUE, ...)
      self$d <- dev.cur()

      addHandlerClicked(mag_plus, function(h, ...) {
        s0 <- gWidgets2::size(gg)
        gWidgets2::size(gg) <- s0*1.1
      })

      addHandlerClicked(mag_minus, function(h, ...) {
        s0 <- gWidgets2::size(gg)
        gWidgets2::size(gg) <- s0/1.1
      })
    },
    ##' @description Set the graphics device to the device contained in the
    ##' magnifier
    devSet = function() {
      dev.set(self$d)
    }
  )
)
