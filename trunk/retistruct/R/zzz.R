zzz.revision <- function() {
  return(as.integer(gsub("Rev: ", "" ,gsub("\\$", "", "$Rev$"))))
}

.First.lib <- function(lib, pkg) {
  library.dynam("retistruct", pkg, lib)
}
