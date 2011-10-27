plot.flat <- function(x, axt="n", ylim=NULL, ...) {
  UseMethod("plot.flat")
}

plot.flat.default <- function(x, axt="n", ylim=NULL, ...) {
}

plot.polar <- function(r, show.grid=TRUE,
                       grid.col="gray", grid.bg="transparent", 
                       grid.int.minor=15, grid.int.major=45, ...) {
  UseMethod("plot.polar")
}

plot.polar.default <- function(r, show.grid=TRUE,
                               grid.col="gray", grid.bg="transparent", 
                               grid.int.minor=15, grid.int.major=45,
                               flip.horiz=FALSE, ...) {
  plot.new()
}

plot.spherical <- function(r, ...) {
  UseMethod("plot.spherical")
}

plot.spherical.default <- function(r, ...) {
  rgl.clear()
  rgl.bg(color="white")
}

plot.l.vs.L <- function(r) {
  UseMethod("plot.l.vs.L")
}

plot.l.vs.L.default <- function(r) {
}

getDss <- function(x) {
  UseMethod("getDss")
}

getDss.default <- function(x) {
  return(NULL)
}

getDss.mean <- function(x) {
  UseMethod("getDss.mean")
}

getDss.mean.default <- function(x) {
  return(NULL)
}

getSss <- function(x) {
  UseMethod("getSss")
}

getSss.default <- function(x) {
  return(NULL)
}

getTss <- function(x) {
  UseMethod("getTss")
}

getTss.default <- function(x) {
  return(NULL)
}
