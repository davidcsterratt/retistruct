.First.lib <- function(lib, pkg) {
  library.dynam("retistruct", pkg, lib)
  options(TF.col="darkcyan")
  options(TB.col="darkcyan")
  options(V.col="darkcyan")
  options(stitch.col="cyan")
  options(V.stitch.col="darkcyan")
}
