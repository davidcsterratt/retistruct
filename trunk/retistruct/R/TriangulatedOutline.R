plot.flat.triangulatedOutline <- function(t, axt="n", ylim=NULL, ...) {
  NextMethod()
  args <- list(...)
  plot.mesh <- is.null(args$mesh) || args$mesh

  if (plot.mesh) 
    with(t, trimesh(T, P, col="grey", add=TRUE))  
}
