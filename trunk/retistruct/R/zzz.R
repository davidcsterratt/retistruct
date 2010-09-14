zzz.revision <- function() {
  return(as.integer(gsub("Rev: ", "" ,gsub("\\$", "", "$Rev: 259$"))))
}

.First.lib <- function(lib, pkg) {
  library.dynam("retistruct", pkg, lib)
}
